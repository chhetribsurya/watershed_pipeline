#!/usr/bin/env Rscript

# Generate rare variant EXPRESSION OUTLIER FILE for later PVALUE conversion

library(data.table)
library(dplyr)

#rv.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt'
#exp.outlier.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputfiles/kgpexV3.GLOBAL.outlier.controls.v3ciseQTLs_globalOutliersRemoved.txt'

# save gene-level annotation file to dir:
#genelevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
#dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# save variant-level annotation file to dir:
#variantlevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
#dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)

rv.file <- args[1]
exp.outlier.file <- args[2]
outfile1 <- args[3]
outfile2 <- args[4]
#pop <- args[7]

# Create directories
#dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rv = fread(rv.file)
exp.outlier = fread(exp.outlier.file)

# Select genes that have at least one "outlier" status
outlier_genes.collapse <- exp.outlier %>% 
  group_by(Gene) %>% 
  filter(any(Status == "outlier")) %>%
  ungroup() %>%
  rename(GeneName = Gene, SubjectID = Individual) %>%
  select(SubjectID, GeneName, Zscore)


# Unique genes list with outlier status
outlier_uniq_geneslist <- exp.outlier %>% filter(Status == "outlier") %>% distinct(Gene)
outlier_genecount <-  nrow(outlier_uniq_geneslist)
cat("\n")
flush.console()
print(paste("Total unique genes with at least one outlier inidividual status:", outlier_genecount))
flush.console()
cat("\n")
flush.console()

# Match rare variants with expression outliers
exp.outlier_rename = exp.outlier %>% rename(GeneName = Gene, SubjectID = Individual)
rv_rename = rv %>% rename(GeneName = Gene, SubjectID = Ind) 

outlier.collapse_init = inner_join(rv_rename, outlier_genes.collapse, by = c("GeneName", "SubjectID"))
outlier.collapse_init$variantID <- paste(outlier.collapse_init$GeneName, paste(outlier.collapse_init$Chrom, outlier.collapse_init$End, outlier.collapse_init$Ref, outlier.collapse_init$Alt, sep = "_"), sep=":")
outlier.collapse_select <- outlier.collapse_init %>% select(SubjectID, GeneName, variantID, Zscore)
outlier_collapse_genecount <-  nrow(outlier.collapse_select)
cat("\n")
flush.console()
print(paste0("Total unique Indiv-gene-rvs with at least one outlier inidividual status for gene:", outlier_collapse_genecount))
flush.console()
cat("\n")
flush.console()

# Like gene level annotation
uniq_outlier_indiv_genes <- outlier.collapse_select %>% distinct(GeneName,SubjectID, .keep_all=T) %>% select(SubjectID,GeneName,Zscore)
uniq_indiv_genecount <- outlier.collapse_select %>% distinct(GeneName,SubjectID) %>% nrow
cat("\n")
flush.console()
print(paste0("Total unique Indiv-gene pairs with at least one rarevariant/outlier inidividual status:", uniq_indiv_genecount))
flush.console()
cat("\n")
flush.console()

#outlier_outfile = paste0(tools::file_path_sans_ext(exp.outlier.file), '.uncollapsed.rvpair.filtgenes.tsv')
write.table(outlier.collapse_select, outfile1, quote = F, sep = '\t', row.names = F, col.names = T)
write.table(uniq_outlier_indiv_genes, outfile2, quote = F, sep = '\t', row.names = F, col.names = T)
print(paste0("Unique Indiv-gene pair with rare vars file saved to:", outfile1))
print(paste0("Unique Indiv-gene pair file saved to:", outfile2))

# copy variant-level annotation file to target annotation dir:
#file.copy(outfile1, output_dir, overwrite = TRUE)


