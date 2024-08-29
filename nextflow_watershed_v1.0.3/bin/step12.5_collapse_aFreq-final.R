#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

rv.file <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]
genelevel_dir <- args[4]
variantlevel_dir <- args[5]


# Create directories
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# Load rarevar file
rv = fread(rv.file)

# Rename columns
rv = rv %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level transformation
rv_init = rv %>% group_by(SubjectID, GeneName) %>%
  summarise(af=min(AF), num_rare_variants=n()) %>%
  select(c("SubjectID", "GeneName", "af", "num_rare_variants"))

#outfile = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.kgpex_annoFreq.collapse.tsv'
write.table(rv_init, outfile1, quote = F, sep = '\t', row.names = F, col.names = T)

# copy gene-level annotation file to target annotation dir:
file.copy(outfile1, genelevel_dir, overwrite=TRUE)

# Variant Level Annotation
rv_filt <- rv %>% rename(af=AF)
rv_filt$variantID <- paste(rv_filt$GeneName, paste(rv_filt$Chrom, rv_filt$End, rv_filt$Ref, rv_filt$Alt, sep = "_"), sep=":")
rv_filt$num_rare_variants <- 1
rv_select <- rv_filt %>% select(SubjectID, GeneName, variantID, af, num_rare_variants)

#rv.outfile = paste0(tools::file_path_sans_ext(rv.file), '.kgpex_annoFreq.uncollapsed.rvpair.tsv')
write.table(rv_select, outfile2, quote = F, sep = '\t', row.names = F, col.names = T)

# copy variant-level annotation file to target annotation dir:
file.copy(outfile2, variantlevel_dir, overwrite=TRUE)



