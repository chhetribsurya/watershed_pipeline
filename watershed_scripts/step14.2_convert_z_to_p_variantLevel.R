#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(stringr)
library(optparse)

# Read command line arguments
option_list = list(make_option(c('--ZFILE'), type = 'character', default = NULL, help = 'File with Z-scores'), 
                   make_option(c('--ZCOL'), type = 'character', default = "MedZ", help = 'Name of column holding zscores'),
                   make_option(c('--OUTFILE'), type = 'character', default = NULL, help = 'Name of desired outputfile'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input.file = opt$ZFILE
output.file = opt$OUTFILE

# open input file
input = fread(input.file)

# select relevant columns
input <- input %>%
  mutate(PVAL = sign(Zscore)*2*pnorm(Zscore, lower.tail=FALSE))

#input <- input %>%
#  mutate(PVAL = 2 * pmin(pnorm(Zscore), pnorm(-Zscore)))

input <- input %>% select(SubjectID, GeneName, variantID, PVAL)
fwrite(input,file=output.file, sep="\t")
