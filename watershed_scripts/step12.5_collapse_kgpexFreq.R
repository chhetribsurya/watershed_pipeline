#!/usr/bin/env Rscript

# COLLAPSE KGPEx rare-variant anno Number and AF

## Collapse GTEx annotations to gene-individual pair level
## The annotations and their transformations across multiple rare variants
## - af (min)
## - num_rare_variants (sum)
## TODO
## - donor_ss (max)
## - acceptor_ss (max)
## - ppt_region (max)
## - donor_ss_window (max)
## - acceptor_ss_winow (max)

library(data.table)
library(dplyr)

rv.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt'
rv = fread(rv.file)

# save gene-level annotation file to dir:
genelevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# save variant-level annotation file to dir:
variantlevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)

# Rename columns
rv = rv %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level transformation
rv_init = rv %>% group_by(SubjectID, GeneName) %>%
  summarise(af=min(AF), num_rare_variants=n()) %>%
  select(c("SubjectID", "GeneName", "af", "num_rare_variants"))

outfile = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.kgpex_annoFreq.collapse.tsv'
write.table(rv_init, outfile, quote = F, sep = '\t', row.names = F, col.names = T)

# copy gene-level annotation file to target annotation dir:
file.copy(outfile, genelevel_dir)

# Variant Level Annotation
rv_filt <- rv %>% rename(af=AF)
rv_filt$variantID <- paste(rv_filt$GeneName, paste(rv_filt$Chrom, rv_filt$End, rv_filt$Ref, rv_filt$Alt, sep = "_"), sep=":")
rv_filt$num_rare_variants <- 1
rv_select <- rv_filt %>% select(SubjectID, GeneName, variantID, af, num_rare_variants)

rv.outfile = paste0(tools::file_path_sans_ext(rv.file), '.kgpex_annoFreq.uncollapsed.rvpair.tsv')
write.table(rv_select, rv.outfile, quote = F, sep = '\t', row.names = F, col.names = T)

# copy variant-level annotation file to target annotation dir:
file.copy(rv.outfile, variantlevel_dir)


