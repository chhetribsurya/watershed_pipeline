#!/usr/bin/env Rscript

# COLLAPSE GENCODE FILE
## Collapse gencode annotations to gene-individual pair level
## Take minimum across multiple rare variants

library(data.table)
library(dplyr)

#gencode.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/gencode/gene-GLOBAL-rv.gencode.txt'

# save gene-level annotation file to dir:
#genelevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
#dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# save variant-level annotation file to dir:
#variantlevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
#dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

gencode.file <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]
genelevel_dir <- args[4]
variantlevel_dir <- args[5]
#phylop.file <- args[6]
#pop <- args[7]

# Create directories
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

gencode = fread(gencode.file)

# Gene Level Annotation
gencode.collapse = gencode %>% group_by(Gene, Ind) %>% 
  summarise(distTSS=min(distTSS), distTES=min(distTES)) %>%
  select("Ind", "Gene", "distTSS", "distTES")

colnames(gencode.collapse) = c("SubjectID", "GeneName", "distTSS", "distTES")

#outfile = paste0(tools::file_path_sans_ext(gencode.file), '.collapse.tsv')
write.table(gencode.collapse, outfile1, quote = F, sep = '\t', row.names = F, col.names = T)

# copy gene-level annotation file to target annotation dir:
file.copy(outfile1, genelevel_dir, overwrite=TRUE)

# Variant Level Annotation
gencode_filt <- gencode %>% rename(GeneName = Gene, SubjectID = Ind)
gencode_filt$variantID <- paste(gencode_filt$GeneName, paste(gencode_filt$Chrom, gencode_filt$End, gencode_filt$Ref, gencode_filt$Alt, sep = "_"), sep=":")
#gencode_filt %>% head
gencode_select <- gencode_filt %>% select(SubjectID, GeneName, variantID, distTSS, distTES)

#gencode.outfile = paste0(tools::file_path_sans_ext(gencode.file), '.uncollapsed.rvpair.tsv')
write.table(gencode_select, outfile2, quote = F, sep = '\t', row.names = F, col.names = T)

# copy variant-level annotation file to target annotation dir:
file.copy(outfile2, variantlevel_dir, overwrite=TRUE)

