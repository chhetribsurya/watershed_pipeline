#!/usr/bin/env Rscript
 
# Converts rare variants file to bed input for querying PhyloP 100way conservation scores from UCSC

library(data.table)
library(dplyr)
library(stringr)
library(optparse)

# Read command line arguments
option_list = list(
  make_option(c('--rv_sites'), type = 'character', default = NULL, help = 'Rare variant file'),
  make_option(c('--outfile'), type = 'character', default = NULL, help = 'Output BED file')
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

rv_file = opt$rv_sites
out_file = opt$outfile

if (is.null(rv_file) || is.null(out_file)) {
  stop("Both --RV and --OUT arguments are required")
}

# open rare variants
rv = fread(rv_file)

# select relevant columns
rv = rv %>% select(Chrom, Start, End)

# convert chromosome to number only
chrom_num = str_remove(rv$Chrom, 'chr')
rv$Chrom = as.numeric(chrom_num)

# sort and remove duplicate rows
rv = rv %>% arrange(Chrom, Start) %>% distinct()

# add the 'chr' prefix back to chromosome column
rv = rv %>% mutate(Chrom = paste0('chr', Chrom))

# write to output file
write.table(rv, out_file, quote = F, sep = '\t', row.names = F, col.names = F)
cat("\n\n Rare variant file converted to bed format for phyloP query ...\n\n")
