#!/usr/bin/env Rscript

# Required libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
#library(gridExtra)
#library(scales)

# Function to display usage instructions
usage <- function() {
  cat("\nUsage: Rscript rare_var_enrichment.R [DATA_DIR] [VARIABLE] [ZSCORE_FILE] [RAREVAR_FILE]\n",
      "Arguments:\n",
      "  OUTPUT_DIR   : Path to the data directory (e.g., /path/to/data).\n",
      "  VARIABLE     : Cohort name for example: 'GLOBAL', 'EUR', or 'AFR'.\n",
      "  ZSCORE_FILE  : Path to the Z-score file.\n",
      "  RAREVAR_FILE : Path to the Rare Variant file.\n\n",
      "Example:\n",
      "  Rscript rare_var_enrichment.R /path/to/data GLOBAL zscore.txt rarevar.txt\n",
      sep="")
  quit(status = 1)
}

cat("\n\n
#######################################################
#                                                     #
# Rare variant enrichment analysis (Zscorebased)      #
#                                                     #
#######################################################
\n\n")


enrichment_compute <- function(outdir, variable, zscore_file, rarevar_file) {

  # Setting file paths and output directory
  z_scores_file <- file.path(zscore_file)
  rare_var_pairs_file <- file.path(rarevar_file)

  # Ensure output directory exists
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Load Zscore file
  gene_individual_data <- fread(z_scores_file, sep="\t")

  # Load and read rare variant file
  rare_variant_data <- fread(rare_var_pairs_file, sep="\t") %>%
    rename(Individual=Ind)

  # Filter genes with at least one outlier individual
  genes_with_outliers <- unique(gene_individual_data[gene_individual_data$Status == "outlier",]$Gene)
  gene_individual_data <- gene_individual_data[gene_individual_data$Gene %in% genes_with_outliers, ]

  # Create a set of gene-individual pairs from the rare variants dataframe
  pairs_with_rare_variants <- rare_variant_data %>%
    mutate(gene_ind_key = paste(Gene, Individual, sep = "_")) %>%
    pull(gene_ind_key) %>% 
    unique()

  #print(head(pairs_with_rare_variants))

  compute_metrics <- function(threshold) {
    print(paste("Processing threshold:", threshold))

    # Identify cases (outliers) and controls based on Zscore
    cases <- gene_individual_data[abs(gene_individual_data$Zscore) >= threshold,]
    controls <- gene_individual_data[abs(gene_individual_data$Zscore) < threshold,]
    
    # Identify which individuals have rare variants
    cases_with_variant <- paste(cases$Gene, cases$Individual, sep = "_") %in% pairs_with_rare_variants
    controls_with_variant <- paste(controls$Gene, controls$Individual, sep = "_") %in% pairs_with_rare_variants
    # Contingency table
    a <- sum(cases_with_variant) # Cases with rare variants
    b <- nrow(cases) - a # Cases without rare variants
    c <- sum(controls_with_variant) # Controls with rare variants
    d <- nrow(controls) - c # Controls without rare variants

    # Convert to numeric to prevent integer overflow
    a <- as.numeric(a)
    b <- as.numeric(b)
    c <- as.numeric(c)
    d <- as.numeric(d)

    # Relative Risk
    rr <- (a / (a + b)) / (c / (c + d))

    # Confidence interval for Relative Risk - using standard error approximation
    rr_se <- sqrt((1/a) - (1/(a+b)) + (1/c) - (1/(c+d)))
    rr_lower <- exp(log(rr) - 1.96 * rr_se)
    rr_upper <- exp(log(rr) + 1.96 * rr_se)

    # Odds Ratio
    or_value <- (a * d) / (b * c)

    # Confidence interval for Odds Ratio - using standard error approximation
    or_se <- sqrt(1/a + 1/b + 1/c + 1/d)
    or_lower <- exp(log(or_value) - 1.96 * or_se)
    or_upper <- exp(log(or_value) + 1.96 * or_se)

    # Fisher's Exact Test for significance of RR and OR
    p_value_rr <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
    p_value_or <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value

    # Chi-squared test for significance of OR
    #p_value_or <- chisq.test(matrix(c(a, b, c, d), nrow = 2))$p.value

    data.frame(Threshold = threshold, 
           Relative_Risk = rr, RR_CI_Lower = rr_lower, RR_CI_Upper = rr_upper, 
           Odds_Ratio = or_value, OR_CI_Lower = or_lower, OR_CI_Upper = or_upper, 
           P_Value_RR = p_value_rr, P_Value_OR = p_value_or,
           Outlier_All_Pairs = a + b, Outlier_Rare_Pairs = a)

  }

  # Calculate for each threshold
  thresholds <- seq(1, 5, 0.5)
  results_df <- bind_rows(lapply(thresholds, compute_metrics))

}

# Function to save results to a file
save_results <- function(results_df, outdir, variable) {
  outfile <- file.path(outdir, paste0("enrichment_analysis_", variable, "_results.tsv"))
  write.table(results_df, outfile, quote = FALSE, sep = '\t', row.names = FALSE, col.names = TRUE)
  return(outfile)
}

# Function to plot Odds Ratio vs Threshold
plot_odds_ratio <- function(results_df, outdir, variable) {
  plot <- ggplot(results_df, aes(x = Threshold, y = Odds_Ratio)) + 
    geom_point(color = 'red') + 
    geom_text(aes(label = Outlier_Rare_Pairs), hjust = 0.7, vjust = -4, size=4) +
    geom_errorbar(aes(ymin = OR_CI_Lower, ymax = OR_CI_Upper), width = 0.1, color = 'red') +
    labs(title = "Odds Ratio vs Threshold", x = "Threshold", y = "Odds Ratio") +
    theme_bw()
  
  plot_file <- file.path(outdir, paste0("odds_ratio_plot_", variable, ".pdf"))
  ggsave(plot_file, plot = plot, width = 8, height = 8)
}

# Function to plot Relative Risk vs Threshold
plot_relative_risk <- function(results_df, outdir, variable) {
  plot <- ggplot(results_df, aes(x = Threshold, y = Relative_Risk)) + 
    geom_point(color = 'red') +
    geom_text(aes(label = Outlier_Rare_Pairs), hjust = 0.7, vjust = -4, size=4) +
    geom_errorbar(aes(ymin = RR_CI_Lower, ymax = RR_CI_Upper), width = 0.1, color = 'red') +
    labs(title = "Relative Risk vs Threshold", x = "Threshold", y = "Relative Risk") +
    theme_bw()

  plot_file <- file.path(outdir, paste0("relative_risk_plot_", variable, ".pdf"))
  ggsave(plot_file, plot = plot, width = 8, height = 8)
}

# Function to create and save a combined plot with OR as a secondary axis
create_combined_plot <- function(results_df, outdir, variable) {
  # Creating a combined plot
  combined_plot <- ggplot(data = results_df, aes(x = Threshold)) +
    geom_point(aes(y = Relative_Risk, color = "Relative Risk")) +
    geom_errorbar(aes(ymin = RR_CI_Lower, ymax = RR_CI_Upper, color = "Relative Risk"), width = 0.1) +
    geom_point(aes(y = Odds_Ratio, color = "Odds Ratio")) +
    geom_errorbar(aes(ymin = OR_CI_Lower, ymax = OR_CI_Upper, color = "Odds Ratio"), width = 0.1) +
    scale_color_manual(values = c("Relative Risk" = "blue", "Odds Ratio" = "red")) +
    labs(y = "Value", color = "Measure") +
    theme_minimal()

  # Saving the plot
  combined_plot_file <- file.path(outdir, paste0("combined_plot_", variable, ".pdf"))
  ggsave(combined_plot_file, plot = combined_plot, width = 10, height = 6)

  return(combined_plot)
}

# Main function
main <- function(output_dir, variable, zscore_file, rarevar_file) {

  # Print the command line arguments
  cat("\nRunning script with the following arguments:\n",
      "OUTPUT_DIR: ", output_dir, "\n",
      "VARIABLE: ", variable, "\n",
      "ZSCORE_FILE: ", zscore_file, "\n",
      "RAREVAR_FILE: ", rarevar_file, "\n\n")

  outdir <- file.path(output_dir, paste0("rarevar_enrichment_", variable))
  results_df <- enrichment_compute(outdir, variable, zscore_file, rarevar_file)
  outfile <- save_results(results_df, outdir, variable)

  plot_odds_ratio(results_df, outdir, variable)
  plot_relative_risk(results_df, outdir, variable)
  #create_combined_plot(results_df, outdir, variable)

}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 4) {
  usage()
}

# Run the main function with provided arguments
main(args[1], args[2], args[3], args[4])

