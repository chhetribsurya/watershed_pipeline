# COLLAPSE PHYLOP FILE

#!/usr/bin/env Rscript

## Collapse PhyloP (UCSC conservation scores) to gene-individual pair level
## Take maximum across multiple rare variants

library(data.table)
library(dplyr)

#rv.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt'
#phylop.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/ucsc/gene-GLOBAL-rv.phyloP100way.bed'

# save gene-level annotation file to dir:
#genelevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
#dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# save variant-level annotation file to dir:
#variantlevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
#dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

rv.file <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]
genelevel_dir <- args[4]
variantlevel_dir <- args[5]
phylop.file <- args[6]
#pop <- args[7]

# Create directories
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# Load rarevar file
rv = fread(rv.file)
phylop = fread(phylop.file, col.names=c("Chrom", "Start", "End", "phylop"))

# Match rare variants with their annotations
phylop.collapse_init = left_join(rv, phylop, by = c("Chrom", "Start"))

# Rename columns
phylop.collapse_init = phylop.collapse_init %>% rename(GeneName = Gene, SubjectID = Ind)

# Gene-Ind level Annotation transformation
phylop.collapse = phylop.collapse_init %>% group_by(GeneName, SubjectID) %>%
  summarise(phylop=max(phylop)) %>%
  relocate("SubjectID")

#outfile = paste0(tools::file_path_sans_ext(phylop.file), '.collapse.tsv')
write.table(phylop.collapse, outfile1, quote = F, sep = '\t', row.names = F, col.names = T)

# copy gene-level annotation file to target annotation dir:
file.copy(outfile1, genelevel_dir, overwrite=TRUE)

# Variant Level Annotation
phylop.collapse_init$variantID <- paste(phylop.collapse_init$GeneName, paste(phylop.collapse_init$Chrom, phylop.collapse_init$End.x, phylop.collapse_init$Ref, phylop.collapse_init$Alt, sep = "_"), sep=":")
phylop.collapse_select <- phylop.collapse_init %>%
    select(SubjectID, GeneName, variantID, phylop)

#phylop.outfile = paste0(tools::file_path_sans_ext(phylop.file), '.uncollapsed.rvpair.tsv')
write.table(phylop.collapse_select, outfile2, quote = F, sep = '\t', row.names = F, col.names = T)

# copy variant-level annotation file to target annotation dir:
file.copy(outfile2, variantlevel_dir, overwrite=TRUE)

