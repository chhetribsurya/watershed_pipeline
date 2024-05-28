#!/usr/bin/env Rscript

# COLLAPSE VEP LOFTEE FILE

# load R
#ml r/4.2.0

## Collapse VEP and LOFTEE annotations to gene-individual pair level
## Take maximum across multiple rare variants

# load R
system("ml r/4.2.0")
library(data.table)
library(dplyr)

rv.file = '/scratch16/abattle4/surya/datasets/WatershedAFR/data/rare_variants_gnomad/gene-GLOBAL-rv.txt'
vep.loftee.file = "/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/vep/gene-GLOBAL-rv.vep.loftee.snv.tsv"

# save gene-level annotation file to dir:
genelevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/genelevel_output"
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

# save variant-level annotation file to dir:
variantlevel_dir="/scratch16/abattle4/surya/datasets/WatershedAFR/data/annotation/annotation_inputFinal/variantlevel_output"
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)

# Load rv and vep files
rv = fread(rv.file)
vep.loftee = fread(vep.loftee.file)

# Convert genomic coordinates in vep.loftee to be consistent with rare variants (0-based)
vep.loftee$Pos = vep.loftee$Pos - 1

# Match rare variants with their annotations
vep.loftee.collapse_init = left_join(rv, vep.loftee, by = c("Chrom" = "Chrom", "Start" = "Pos"))

# Rename columns
vep.loftee.collapse_init = vep.loftee.collapse_init %>% rename(GeneName = Gene, SubjectID = Ind)


# Gene-Individual level transformation
vep.loftee.collapse_grouped = vep.loftee.collapse_init %>% group_by(SubjectID, GeneName) %>%
  summarise(`3_prime_UTR_variant`=max(`3_prime_UTR_variant`), `5_prime_UTR_variant`=max(`5_prime_UTR_variant`), 
            TF_binding_site_variant=max(TF_binding_site_variant), downstream_gene_variant=max(downstream_gene_variant),
            intergenic_variant=max(intergenic_variant), intron_variant=max(intron_variant), missense_variant=max(missense_variant),
            non_coding_transcript_exon_variant=max(non_coding_transcript_exon_variant),
            regulatory_region_variant=max(regulatory_region_variant), splice_acceptor_variant=max(splice_acceptor_variant),
            splice_donor_variant=max(splice_donor_variant), splice_region_variant=max(splice_region_variant),
            stop_gained=max(stop_gained), synonymous_variant=max(synonymous_variant),
            upstream_gene_variant=max(upstream_gene_variant), LoF_HC=max(LoF_HC), LoF_LC=max(LoF_LC)) %>%
  relocate("SubjectID")

outfile = paste0(tools::file_path_sans_ext(vep.loftee.file), '.collapse.tsv')
write.table(vep.loftee.collapse_grouped, outfile, quote = F, sep = '\t', row.names = F, col.names = T)

# copy gene-level annotation file to target annotation dir:
file.copy(outfile, genelevel_dir)

######
# Gene-Individual-Variant level transformation i.e, uncollapsed version
vep.loftee.filt <- vep.loftee %>% distinct(Chrom, Pos, .keep_all = TRUE)

vep.loftee.collapse <- left_join(rv, vep.loftee.filt, by = c("Chrom" = "Chrom", "Start" = "Pos"))
vep.loftee.collapse <- vep.loftee.collapse %>% rename(GeneName = Gene, SubjectID = Ind)

vep.loftee.collapse$variantID <- paste(vep.loftee.collapse$GeneName, paste(vep.loftee.collapse$Chrom, vep.loftee.collapse$End, vep.loftee.collapse$Ref.x, vep.loftee.collapse$Alt.x, sep = "_"), sep=":")
vep.loftee.collapse_select <- vep.loftee.collapse %>% select(SubjectID, GeneName, variantID, starts_with('3_prime_UTR_variant'):LoF_LC)

outfile = paste0(tools::file_path_sans_ext(vep.loftee.file), '.uncollapsed.rvpair.tsv')
write.table(vep.loftee.collapse_select, outfile, quote = F, sep = '\t', row.names = F, col.names = T)

# copy variant-level annotation file to target annotation dir:
file.copy(outfile, variantlevel_dir)

# SANITY CHECK 
# Check if all values beyond first four columns are the same for all duplicated rows based on Chrom and Pos
# If the values are same and thus distinct(Chrom, Pos, .keep_all = TRUE) is warranted to remove multi-allelic sites
# For example: vep.loftee.filt <- vep.loftee %>%
#   distinct(Chrom, Pos, .keep_all = TRUE)

setDT(vep.loftee)
vep.loftee_chr22 <- vep.loftee[Chrom == "chr22"]
duplicated_rows <- vep.loftee_chr22[,.SD[.N > 1], by = .(Chrom, Pos)]

result <- duplicated_rows %>%
  group_by(Chrom, Pos) %>%
  mutate(if_all(-c(1:4), ~length(unique(.)) > 1)) %>%
  ungroup() %>%
  filter(if_all(-c(1:4), ~. == TRUE))

# Check if there are any rows in the result
if(nrow(result) > 0) {
  print("There are duplicated rows with differing values beyond the first four columns.")
} else {
  print("All duplicated rows have the same values beyond the first four columns.")
  print("Okay to proceed!")
}

