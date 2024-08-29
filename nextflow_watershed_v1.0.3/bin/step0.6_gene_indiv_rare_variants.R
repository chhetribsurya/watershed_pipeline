#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(data.table)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--rv_sites'), type = 'character', default = NULL, help = 'Bed file with rare variants overlapping 10kb +/- window around the genes'),
                   make_option(c('--popname'), type = 'character', default = NULL, help = 'Population name for output file name'),
                   make_option(c('--outfile'), type = 'character', default = NULL, help = 'Save file with rare variants per gene-individual pair'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

rv_sites_file = opt$rv_sites
popname = opt$popname
outputfile = opt$outfile

# Open file
rv_sites = read.table(rv_sites_file, header = FALSE, check.names = FALSE)

# Rearrange columns
df = rv_sites[,c("V11", "V7", "V8", "V2", "V3", "V4", "V5", "V6")]
colnames(df) = c("Gene","Ind","Chrom","Start","End","Ref","Alt", "AF")

# Sort by genes
df = arrange(df, Gene)

# Save output files
write.table(df, outputfile, sep = '\t', row.names = FALSE, quote = FALSE)
cat(paste0("Output file saved to : ", outputfile, '\n'))
