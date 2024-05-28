#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(data.table)
library(optparse)

# Function to parse command line arguments
parse_command_line_arguments <- function() {
  option_list <- list(
    make_option(c("--rv_sites"), type = "character", default = NULL, help = "bed file with rare variants overlapping 10kb +/- window around the genes"),
    make_option(c("--popname"), type = "character", default = NULL, help = "population of individuals for output file name"),
    make_option(c("--outdir"), type = "character", default = NULL, help = "directory to save list of rare variants per gene-individual pair")
  )
  
  opt_parser <- OptionParser(option_list = option_list)
  parse_args(opt_parser)
}

# Function to read and process rare variant sites
process_rare_variant_sites <- function(rv_sites_file) {
  rv_sites <- read.table(rv_sites_file, header = FALSE, check.names = FALSE)
  
  df <- rv_sites[, c("V11", "V7", "V8", "V2", "V3", "V4", "V5", "V6")]
  colnames(df) <- c("Gene", "Ind", "Chrom", "Start", "End", "Ref", "Alt", "AF")
  
  arrange(df, Gene)
}

# Function to save gene-individual rare variants to file
save_gene_individual_rare_variants <- function(df, outdir, popname) {
  filename <- paste0(outdir, "/gene-", popname, "-rv.txt")
  write.table(df, filename, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(paste0("Saved to ", filename, "\n"))
}

# Main function
main <- function() {
  opt <- parse_command_line_arguments()
  
  rv_sites_file <- opt$rv_sites
  popname <- opt$popname
  outdir <- opt$outdir
  
  df <- process_rare_variant_sites(rv_sites_file)
  save_gene_individual_rare_variants(df, outdir, popname)
}

main()

