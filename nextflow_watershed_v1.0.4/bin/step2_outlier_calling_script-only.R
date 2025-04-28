#!/usr/bin/env Rscript

###################################
#
# CALL Expression OUTLIERS
#
###################################

library(data.table)
library(dplyr)
library(tidyr)

#' Calculate P-values and Assign Status
#'
#' This function calculates p-values from Z-scores, retains the sign, and assigns a status 
#' ("outlier" or "control") based on either Z-score or p-value thresholds.
#'
#' @param df A dataframe in long format with columns: Gene, Individual, Zscore.
#' @param pval_significance_level The overall significance level for p-value calculation.
#' @param z_threshold Threshold for Z-score to determine outliers.
#' @return A modified dataframe with additional columns: PVAL, SignPVAL, and Status.
calculate_pvalues <- function(df, pval_significance_level, z_threshold) {
  #bonferroni_correction <- pval_significance_level / length(unique(df$Individual))
  bonferroni_correction <- pval_significance_level /length(unique(df$Gene))
  
  if (!is.null(z_threshold)) {
    df %>% mutate(PVAL = 2 * pnorm(-abs(Zscore)),
                  SignPVAL = sign(Zscore) * PVAL,
                  Status = ifelse(abs(Zscore) > z_threshold, "outlier", "control"))
  } else {
    df %>% mutate(PVAL = 2 * pnorm(-abs(Zscore)),
                  SignPVAL = sign(Zscore) * PVAL,
                  Status = ifelse(PVAL < bonferroni_correction, "outlier", "control"))
  }
}


#' Remove Global Outliers
#'
#' This function removes individuals identified as global outliers based on the specified criteria.
#' Global outliers are determined by either the proportion or frequency of outliers per individual.
#'
#' @param df A dataframe in long format with columns: Gene, Individual, Zscore, Status.
#' @param global_method Method for identifying global outliers.
#'                     'proportion' to use the proportion of outliers,
#'                     otherwise frequency will be used.
#' @return A dataframe with global outliers removed.
remove_global_outliers <- function(df, global_method) {
  zCount <- df %>%
    group_by(Individual) %>%
    summarize(OutlierCount = sum(Status == "outlier"), 
              Total = n(), 
              PropOut = OutlierCount / Total)

  if (global_method == 'proportion') {
    q3 <- quantile(zCount$PropOut, type = 7)[4]
    q1 <- quantile(zCount$PropOut, type = 7)[2]
  } else {
    q3 <- quantile(zCount$OutlierCount, type = 7)[4]
    q1 <- quantile(zCount$OutlierCount, type = 7)[2]
  }

  qthres <- q3 + 1.5 * (q3 - q1)
  indsRemove <- if (global_method == 'proportion') {
                  unique(zCount %>% filter(PropOut > qthres) %>% pull(Individual))
                } else {
                  unique(zCount %>% filter(OutlierCount > qthres) %>% pull(Individual))
                }

  df_filtered <- df %>% filter(!(Individual %in% indsRemove))
  return(df_filtered)
}

calculate_proportion_zscore <- function(data, z_threshold) {
  proportion_above_threshold <- sum(abs(data$Zscore) > z_threshold) / nrow(data)

  cat("\nProportion of Z-scores above the threshold of", z_threshold, ":", 
      proportion_above_threshold, "\n")

  return(proportion_above_threshold)
}


calculate_proportion_pvalue <- function(data, pval_significance_level) {
  #number_of_individuals <- length(unique(data$Individual))
  #bonferroni_corrected_threshold <- pval_significance_level / number_of_individuals
  bonferroni_corrected_threshold <- pval_significance_level /length(unique(data$Gene))

  proportion_below_corrected_threshold <- sum(data$PVAL < bonferroni_corrected_threshold) / nrow(data)

  cat("\nProportion of entries below the Bonferroni-corrected PVAL threshold of", 
      bonferroni_corrected_threshold, ":", 
      proportion_below_corrected_threshold, "\n")

  return(proportion_below_corrected_threshold)
}


#' Process Expression Outlier Data
#'
#' Processes eOutlier signal data by converting Z-scores to p-values, assigning outlier or control status,
#' and removing global outliers. The function allows using either Z-score threshold or 
#' p-value significance level for determining the status.
#'
#' @param zscore_data A zscore matrix with eOutlier data where the first column is 'Gene' and the 
#'             subsequent columns are individuals with their Z-scores.
#' @param z_threshold (Optional) A numeric value setting the Z-score threshold to determine outliers.
#'                    If NULL, p-value significance level is used. Default is NULL.
#' @param pval_significance_level The overall desired significance level for p-value thresholding.
#' @param global_method A character string specifying the method to remove global outliers.
#' @return A dataframe with columns: Gene, Individual, Status, Zscore, PVAL, and SignPVAL,
#'         with global outliers removed.
process_outlier_data <- function(zscore_data, z_threshold = NULL, pval_significance_level, global_method) {

  # Read the data
  data <- read.table(zscore_data, header = TRUE, stringsAsFactors = FALSE, sep = "\t")

  long_format <- data %>% 
    pivot_longer(cols = -Id, names_to = "Individual", values_to = "Zscore") %>%
    rename(Gene = Id)

  # Calculate p-values and assign status
  long_format <- calculate_pvalues(long_format, pval_significance_level, z_threshold)

  # Determine the method for status assignment and calculate top fraction
  if (!is.null(z_threshold)) {
    calculate_proportion_zscore(long_format, z_threshold)
  } else {
    calculate_proportion_pvalue(long_format, pval_significance_level)
  }

  # Remove global outliers
  long_format <- remove_global_outliers(long_format, global_method)

  # Rearrange columns
  final_data <- long_format %>% 
    select(Gene, Individual, Status, Zscore, PVAL, SignPVAL)

  return(final_data)
}


#' Process Outlier Genes for filters
#'
#' This function processes outlier genes by joining them with rare variants data
#' and selecting distinct gene-individual pairs with their Z-scores.
#'
#' @param final_data A dataframe with columns: Gene, Individual, Status, Zscore, PVAL, SignPVAL.
#' @param datadir Path to the data directory.
#' @param rv_filename Relative path to the rare variants file within the data directory.
#' @return A dataframe of unique outlier individuals and genes with their Z-scores.
process_outlier_genes <- function(final_data, rv_file) {

  # Select genes that have at least one "outlier" status
  outlier_genes_collapse <- final_data %>%
    group_by(Gene) %>%
    filter(any(Status == "outlier")) %>%
    ungroup() %>%
    rename(GeneName = Gene, SubjectID = Individual) %>%
    select(SubjectID, GeneName,  Status, Zscore, PVAL, SignPVAL)

  # Read rare variants data
  rv <- fread(rv_file)
  rv_rename <- rv %>% rename(GeneName = Gene, SubjectID = Ind)

  # Join with rare variants data
  outlier_collapse_init <- inner_join(rv_rename, outlier_genes_collapse, by = c("GeneName", "SubjectID"))
  outlier_collapse_init$variantID <- paste(outlier_collapse_init$GeneName, paste(outlier_collapse_init$Chrom, outlier_collapse_init$End, outlier_collapse_init$Ref, outlier_collapse_init$Alt, sep = "_"), sep=":")
  outlier_collapse_select <- outlier_collapse_init %>% select(SubjectID, GeneName, variantID, Zscore, PVAL, SignPVAL)

  # Select distinct gene-individual pairs
  uniq_outlier_indiv_genes <- outlier_collapse_init %>%
    distinct(GeneName, SubjectID, .keep_all = TRUE) %>%
    select(SubjectID, GeneName, Zscore, PVAL, SignPVAL)


  # Return both dataframes as elements of a list
  return(list(outlier_collapse_select = outlier_collapse_select, 
              uniq_outlier_indiv_genes = uniq_outlier_indiv_genes))
}


# Usage
# Main function to run the script
#main <- function(expr_file, rv_file, output_dir, pop) {
main <- function(expr_file, output_dir, pop) {
    #datadir="/scratch16/abattle4/surya/datasets/WatershedAFR/data"
    #zscore_data <- file.path(datadir, "/data_prep/PEER_EUR/LCL_Blood.peer.v3ciseQTLs.ztrans.txt")
    #rv_file = file.path(datadir, "rare_variants_gnomad/EUR/gene-EUR-rv.txt")
    
    # Set output dir
    #outdir <- file.path(datadir, "rv_expoutlier_refbase/EUR")

    zscore_data <- expr_file
    #rv_file <- rv_file
    outdir <- output_dir

    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    # Function parameters
    pval_significance_level <- 0.05  # Desired overall significance level (e.g., 0.05)
    global_method <- 'proportion'  # Method for global outlier removal
    zthreshhold <- 3

    cat("\n
    ###########################
    #                         #
    # CALLING  EXP. OUTLIERS  #
    #                         #
    ###########################
    \n")

    # Based on Pvalue threshold outlier calls
    finaloutfile_pval <- file.path(outdir, paste0("Cohort.", pop, ".Pvalthresbased.globalOutliersRemoved.txt"))
    final_data_pval <- process_outlier_data(zscore_data, pval_significance_level = pval_significance_level, global_method = global_method)

    write.table(final_data_pval, finaloutfile_pval, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nData saved to", finaloutfile_pval, "\n")

    # Based on Zscore threshold outlier calls
    finaloutfile_zscore <- file.path(outdir, paste0("Cohort.", pop, ".Zscorethresbased.globalOutliersRemoved.txt"))
    final_data_zscore <- process_outlier_data(zscore_data, pval_significance_level = pval_significance_level, z_threshold=zthreshhold, global_method = global_method)

    write.table(final_data_zscore, finaloutfile_zscore, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("\nData saved to", finaloutfile_zscore, "\n")
    cat("\nGLOBAL Outlier Calls Completed ...\n")

}

    # Generate PVAL-based reference-base expOutliers with rare variants presence
    #result_pval <- process_outlier_genes(final_data_pval, rv_file)
    #gene_indiv_unique_df_pval <- result_pval$uniq_outlier_indiv_genes
    #gene_indiv_withrvs_df_pval <- result_pval$outlier_collapse_select

    ## Write Pval threshold based files
    #outfile_pval1 = file.path(outdir, paste0("gene-", pop, "-rv.expoutlier.Pvalthresbased.reference.tsv"))
    #outfile_pval2 = file.path(outdir, paste0("gene-", pop, "-rv.expoutlier.withrvIds.Pvalthresbased.reference.tsv"))
    #write.table(gene_indiv_unique_df_pval, outfile_pval1, sep = "\t", row.names = FALSE, quote = FALSE)
    #write.table(gene_indiv_withrvs_df_pval, outfile_pval2, sep = "\t", row.names = FALSE, quote = FALSE)
    #cat("\nMain Data saved to", outfile_pval1, "\n")

    ## Generate ZSCORE-based reference-base expOutliers with rare variants presence
    #result_zscore <- process_outlier_genes(final_data_zscore, rv_file)
    #gene_indiv_unique_df_zscore <- result_zscore$uniq_outlier_indiv_genes
    #gene_indiv_withrvs_df_zscore <- result_zscore$outlier_collapse_select

    ## Write Zthreshold based files
    #outfile_zscore1 = file.path(outdir, paste0("gene-", pop, "-rv.expoutlier.Zscorethresbased.reference.tsv"))
    #outfile_zscore2 = file.path(outdir, paste0("gene-", pop, "-rv.expoutlier.withrvIds.Zscorethresbased.reference.tsv"))
    #write.table(gene_indiv_unique_df_zscore, outfile_zscore1, sep = "\t", row.names = FALSE, quote = FALSE)
    #write.table(gene_indiv_withrvs_df_zscore, outfile_zscore2, sep = "\t", row.names = FALSE, quote = FALSE)
    #cat("\nMain Data saved to", outfile_zscore1, "\n")
    #cat("\nGLOBAL Outlier Calls Processing Completed for N2pairsInput...\n")

#}

# Function usage
usage <- function() {
  cat("\n
Usage: Rscript outlier_calling_script.R [EXPR_FILE] [OUTPUT_DIR] [POP]\n\n
Arguments:\n
  EXPR_FILE  : Expression file after PEER runs.\n
  OUTPUT_DIR : Path to the output directory.\n
  POP        : Population group, for example: 'GLOBAL', 'EUR', 'AFR'.\n\n
  
Example:\n
  Rscript step2_outlier_calling_script-final.R expOutlier_file /path/to/output/dir GLOBAL\n")
  quit(status = 1)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 3) {
  usage()
}

# Run the main function with provided arguments
main(args[1], args[2], args[3])



