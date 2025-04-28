#!/usr/bin/env R

library(data.table)
library(dplyr)
library(ggplot2)
#library(PRROC)

################################
#                              #
# Watershed Evaluation Script  #
#                              #
################################

# Function usage
usage <- function() {
  cat("\n
Usage: Rscript watershed_evaluation_script.R [EVAL_RDS] [PRED_RDS] [ANNOT_FILE] [OUTPUT_DIR] [POP]\n\n
Arguments:\n
  EVAL_RDS   : Path to the evaluation RDS file.\n
  PRED_RDS   : Path to the prediction RDS file.\n
  ANNOT_FILE : Path to the merged annotation file.\n
  OUTPUT_DIR : Directory to save the plots.\n
  POP        : Population group, e.g., 'GLOBAL', 'EUR', 'AFR'.\n\n
Example:\n
  Rscript watershed_evaluation_script.R eval.rds pred.rds annot.tsv /path/to/output/dir GLOBAL\n")
  quit(status = 1)
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 5) {
  usage()
}

# Function to generate Precision-Recall Curve and save to PDF
#' @description This function generates a Precision-Recall (PR) curve based on watershed and GAM recall and precision data.
#' @param eval_rds Path to the evaluation RDS file.
#' @param output_dir Directory to save the plot.
#' @param pop Population group (e.g., "GLOBAL", "EUR", "AFR").
generate_pr_curve <- function(eval_rds, output_dir, pop) {
  
  # Load Evaluation RDS object
  evaluation_object <- readRDS(eval_rds)
  
  # Watershed PR data for eOutlier
  watershed_recall <- evaluation_object$auc[[1]]$evaROC$watershed_recall
  watershed_precision <- evaluation_object$auc[[1]]$evaROC$watershed_precision
  
  # GAM PR data for eOutlier
  gam_recall <- evaluation_object$auc[[1]]$evaROC$GAM_recall
  gam_precision <- evaluation_object$auc[[1]]$evaROC$GAM_precision
  
  # Create a data frame to hold the data for ggplot2
  watershed_data <- data.frame(recall = watershed_recall, precision = watershed_precision, model = "Watershed-SNV")
  gam_data <- data.frame(recall = gam_recall, precision = gam_precision, model = "GAM-SNV")
  combined_data <- rbind(watershed_data, gam_data)

  # Save the precision-recall data to a CSV file
  pr_data_outfile <- file.path(output_dir, paste0(pop, "_precision_recall_data.csv"))
  write.csv(combined_data, file = pr_data_outfile, row.names = FALSE)
  
  # Function to compute AUPR using the trapezoidal rule
  compute_aupr <- function(recall, precision) {
      # Ensure recall is sorted in ascending order for trapezoidal rule
      sorted_indices <- order(recall)
      recall <- recall[sorted_indices]
      precision <- precision[sorted_indices]
      
      # Compute AUPR using the trapezoidal rule
      aupr <- sum((recall[-1] - recall[-length(recall)]) * ((precision[-1] + precision[-length(precision)]) / 2))
      
      return(aupr)
  }

  # Compute AUPR for Watershed and GAM models
  #pr_watershed <- pr.curve(scores.class0 = watershed_precision, weights.class0 = watershed_recall, curve = TRUE)
  #pr_gam <- pr.curve(scores.class0 = gam_precision, weights.class0 = gam_recall, curve = TRUE)
  
  # Get AUPR values
  #aupr_watershed <- pr_watershed$auc.integral
  #aupr_gam <- pr_gam$auc.integral

  # Compute AUPR for Watershed and GAM models
  aupr_watershed <- compute_aupr(watershed_recall, watershed_precision)
  aupr_gam <- compute_aupr(gam_recall, gam_precision)

  # Generate the ggplot object
  p <- ggplot(combined_data, aes(x = recall, y = precision, color = model, linetype = model)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Watershed-SNV" = "blue", "GAM-SNV" = "red")) +
    scale_linetype_manual(values = c("Watershed-SNV" = "solid", "GAM-SNV" = "dotted")) +
    labs(title = paste("Precision-Recall Curve for", pop), x = "Recall", y = "Precision") +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 2)))
  
  # Save the ggplot object to a PDF
  pr_curve_outfile <- file.path(output_dir, paste0(pop, "_pr_curve_expOutlier.pdf"))
  ggsave(filename = pr_curve_outfile, plot = p, device = "pdf", width = 8, height = 6)

  # Generate the ggplot object
  q <- ggplot(combined_data, aes(x = recall, y = precision, color = model, linetype = model)) +
    geom_line(linewidth = 1) +  # Changed from size to linewidth
    scale_color_manual(values = c("Watershed-SNV" = "blue", "GAM-SNV" = "red")) +
    scale_linetype_manual(values = c("Watershed-SNV" = "solid", "GAM-SNV" = "dotted")) +
    labs(
      title = paste("Precision-Recall Curve for", pop), 
      x = "Recall", y = "Precision",
      subtitle = paste0("AUPR (Watershed-SNV): ", round(aupr_watershed, 3), 
                        " | AUPR (GAM-SNV): ", round(aupr_gam, 3))
    ) +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    guides(color = guide_legend(override.aes = list(linewidth = 2)))
  
  # Save the ggplot object to a PDF
  pr_curve_outfile <- file.path(output_dir, paste0(pop, "_pr_curve_expOutlier_withAUC.pdf"))
  ggsave(filename = pr_curve_outfile, plot = q, device = "pdf", width = 8, height = 6)

}


# Function to generate annotation weights barplot and save to PDF
#' @description This function generates a barplot of annotation weights from the prediction RDS object and saves it as a PDF.
#' @param pred_rds Path to the prediction RDS file.
#' @param annot_file Path to the merged annotation file.
#' @param output_dir Directory to save the plot.
#' @param pop Population group (e.g., "GLOBAL", "EUR", "AFR").
generate_annotation_weights_barplot <- function(pred_rds, annot_file, output_dir, pop) {
  
  # Load prediction RDS object
  prediction_object <- readRDS(pred_rds)
  
  # Extract annotation weights and create a data frame
  annot_weights <- prediction_object$model_params$theta[, 1]
  df_annotation <- fread(annot_file)
  
  # Extract annotation column names, excluding first 2 (indiv, genename) and last 2 (pval, n2pair)
  annot_ids <- colnames(df_annotation)[3:(ncol(df_annotation) - 2)]
  
  # Create a data frame for annotation weights
  df_result <- data.frame(annot_id = annot_ids, weight = annot_weights)
  df_prediction_result <- df_result %>% arrange(desc(weight))
  
  # Adjust factor levels for annotation IDs
  df_prediction_result$annot_id <- factor(df_prediction_result$annot_id,
                                          levels = df_prediction_result$annot_id[order(df_prediction_result$weight)])
  
  # Plot annotation weights
  p <- ggplot(df_prediction_result, aes(x = annot_id, y = weight)) +
    geom_bar(stat = "identity", fill = "blue", color = "black") +
    coord_flip() +  # Flip coordinates for horizontal bars
    labs(title = paste("Annotation Weights for", pop),
         y = "Weight", x = "Annotation") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 6))  # Reduce font size for y-axis labels
  
  # Save the ggplot object to a PDF
  annot_weights_outfile <- file.path(output_dir, paste0(pop, "_annotation_weights_barplot.pdf"))
  ggsave(filename = annot_weights_outfile, plot = p, device = "pdf", width = 8, height = 6)
}


# Main function to run the watershed evaluation and plot generation
main <- function(eval_rds, pred_rds, annot_file, output_dir, pop) {
  
  # Function to generate Precision-Recall Curve
  generate_pr_curve(eval_rds, output_dir, pop)
  
  # Function to generate Annotation Weights Barplot
  generate_annotation_weights_barplot(pred_rds, annot_file, output_dir, pop)
}

# Run the main function
main(args[1], args[2], args[3], args[4], args[5])

