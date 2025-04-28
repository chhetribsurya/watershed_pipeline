#!/usr/bin/env Rscript

# Rscript to process TPM and read count matrices so that columns match covariate files.
# Prepare matrices for input to PEER.

# Load required packages
library(data.table)
library(stringr)
library(optparse)

#' Transform and Normalize Expression Data
#'
#' This function performs sanity checks, subsets and filters TPM and read count matrices
#' based on covariate information, applies log transformation and z-score normalization.
#' It is designed to preprocess gene expression data for further statistical analysis.
#'
#' @param pop A character string specifying the population/cohort type being analyzed.
#' @param tissue A character string specifying the tissue type being analyzed.
#' @param tpm_matrix A data.table object containing the TPM values for genes across samples.
#' @param read_count_matrix A data.table object containing read counts for genes across samples.
#' @param covariates A data.frame containing covariate information with at least one column 'SUBJID' that matches the column names of `tpm_matrix`.
#' @param min_reads Integer specifying the minimum read count threshold for filtering genes.
#' @param min_tpm Numeric specifying the minimum TPM threshold for filtering genes.
#' @return Writes the normalized and transformed data to a file in the current directory and returns NULL.
transform_and_normalize_data = function(pop, tissue, tpm_matrix, read_count_matrix, covariates, min_reads = 6, min_tpm = 0.1) {

    # Sanity checks
    stopifnot(sum(colnames(tpm_matrix) != colnames(read_count_matrix)) == 0)
    stopifnot(sum(tpm_matrix$Gene != read_count_matrix$Gene) == 0)
    print("Sanity check completed for equal matrix of TPM and read count file ...")
    #cat("\n")

    covariates_subset = covariates$SUBJID[covariates$SUBJID %in% colnames(tpm_matrix)]
    genes = tpm_matrix$Gene
    tpm_matrix = tpm_matrix[, covariates_subset, with = FALSE]
    read_count_matrix = read_count_matrix[, covariates_subset, with = FALSE]
    threshold = round(0.2 * ncol(tpm_matrix))
    minimum_present = round(0.05 * ncol(tpm_matrix))
    if (minimum_present == 0) {
        minimum_present = 1
    }
    valid_indices = rowSums(tpm_matrix > min_tpm & read_count_matrix > min_reads) >= threshold
    present_indices = rowSums(tpm_matrix > 0) >= minimum_present
    valid_indices = valid_indices[which(valid_indices %in% present_indices)]
    tpm_matrix = tpm_matrix[valid_indices, ]
    normalized_matrix = scale(t(log2(tpm_matrix + 2)))
    colnames(normalized_matrix) = genes[valid_indices]
    write.table(normalized_matrix, paste0(tissue, ".", pop,  '.log2.ztrans.txt'), quote = FALSE, sep = '\t', row.names = TRUE, col.names = TRUE)
    return()
}


# Command line arguments
option_list = list(
  make_option(c('--COV'), type = 'character', default = NULL, 
              help = 'Path to the covariate file [REQUIRED]', metavar='file'),
  make_option(c('--TPM_FILE'), type = 'character', default = NULL, 
              help = 'Filename for the TPM matrix [REQUIRED]', metavar='file'),
  make_option(c('--READ_FILE'), type = 'character', default = NULL, 
              help = 'Filename for the read count matrix [REQUIRED]', metavar='file'),
  make_option(c('--min_reads'), type = 'integer', default = 6,
              help = 'Minimum read count filter', metavar='integer'),
  make_option(c('--min_tpm'), type = 'numeric', default = 0.1,
              help = 'Minimum TPM filter', metavar='numeric'),
  make_option(c('--pop'), type = 'character', default = 'GLOBAL',
              help = 'Population to be analyzed [REQUIRED]', metavar='numeric'),
  make_option(c('--tissue'), type = 'character', default = 'Blood',
              help = 'Tissue type to be analyzed [REQUIRED]', metavar='string')
)


# Parse arguments
opt_parser = OptionParser(usage = 'Usage: %prog [OPTIONS] <arguments>',
                          option_list = option_list)
opt = parse_args(opt_parser)

# Check if the required arguments are specified
if (is.null(opt$COV) | is.null(opt$TPM_FILE) | is.null(opt$READ_FILE)) {
    stop("The --COV, --TPM_FILE, --READ_FILE, --pop and --tissue options must be specified.\nRun with --help for usage.")
}


# Read in covariates and matrices
covariates = read.table(opt$COV, header = TRUE)
tpm_matrix = fread(opt$TPM_FILE)
read_count_matrix = fread(opt$READ_FILE)
min_reads_thresh = opt$min_reads
min_tpm_thresh = opt$min_tpm
tissue_type = opt$tissue
pop_type = opt$pop

# Call transformation function
transform_and_normalize_data(pop = pop_type, tissue = tissue_type, tpm_matrix, read_count_matrix, covariates, min_reads = min_reads_thresh, min_tpm = min_tpm_thresh)

