#!/usr/bin/env Rscript

# COLLAPSE CADD FILE -- Same as above but Chromosome wise
## Collapse CADD annotations to gene-individual pair level
## Take maximum across multiple rare variants

#system("ml r/4.2.0")
library(data.table)
library(dplyr)
library(tidyr)
library(R.utils)

# To read gz and bz2 files directly, 
# fread() requires 'R.utils' package
# Please install 'R.utils' using 'install.packages('R.utils')'.
# install.packages('R.utils')

args <- commandArgs(trailingOnly = TRUE)

rv.file <- args[1]
outfile1 <- args[2]
outfile2 <- args[3]
genelevel_dir <- args[4]
variantlevel_dir <- args[5]
cadd_base_dir <- args[6]  # base directory for cadd files
pop <- args[7]

# Create directories
dir.create(variantlevel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(genelevel_dir, recursive = TRUE, showWarnings = FALSE)

rv = fread(rv.file)
#header = TRUE

for (chrom in 1:22){
  #chrom=20
  cat("\n\n")
  print(paste0("Processing: chr", chrom))
  
  cadd.file = paste0(cadd_base_dir, paste0('/gene-', pop, '-rv.CADD.chr'), chrom, '.tsv.gz')
  
  cadd = fread(cadd.file) %>%
    rename("Chrom" = "#Chrom") %>%
    mutate(Chrom = paste0("chr",Chrom)) %>%
    rename(GeneSymbol = GeneName) # Change GeneName to GeneSymbol
    #select(-GeneName) #These are Gene symbol and could cause conflict error with rv.file GeneName

  # duplicated_rows <- cadd %>%
  # filter(duplicated(Chrom, Pos) | duplicated(Chrom, Pos, fromLast = TRUE))

  cadd <- cadd %>%
        mutate(SIFTcat = ifelse(is.na(SIFTcat), "missing", as.character(SIFTcat)),
        PolyPhenCat = ifelse(is.na(PolyPhenCat), "missing", as.character(PolyPhenCat)))

  cadd <- cadd %>%
        # Create binary columns for 'SIFTcat'
        pivot_wider(names_from = SIFTcat, values_from = SIFTcat, 
                    names_prefix = "SIFTcat_") %>%
        mutate(across(starts_with("SIFTcat_"), ~ifelse(is.na(.), 0, 1))) %>%
        # Create binary columns for 'PolyPhenCat'
        pivot_wider(names_from = PolyPhenCat, values_from = PolyPhenCat, 
                    names_prefix = "PolyPhenCat_") %>%
        mutate(across(starts_with("PolyPhenCat_"), ~ifelse(is.na(.), 0, 1)))

   #if (!"PolyPhenCat_unknown" %in% names(cadd)) {
   #   cadd$PolyPhenCat_unknown <- 0  # Or another default value as per your requirements
   #}

   # Check if columns exist and set them to 0 if they don't
   columns_to_check <- c("SIFTcat_deleterious", "SIFTcat_tolerated", 
                          "PolyPhenCat_benign", "PolyPhenCat_possibly_damaging", 
                          "PolyPhenCat_unknown")
   for (col in columns_to_check) {
      if (!col %in% names(cadd)) {
        print(paste0(col, " not present and adding a column"))
        cadd[[col]] <- 0
        }
    }

   cadd <- cadd %>%
        rename(PolyPhenCat_NA = PolyPhenCat_missing,
        SIFTcat_NA = SIFTcat_missing)


   # Convert genomic coordinates in CADD annotations to be consistent with rare variants (0-based)
   cadd$Pos = cadd$Pos - 1

   # Match rare variants with their annotations and rename columns
   cadd.collapse_init = filter(rv, Chrom == paste0("chr",chrom)) %>%
    left_join(., cadd, by = c("Chrom" = "Chrom", "Start" = "Pos", "Ref"= "Ref", "Alt"="Alt")) %>%
    rename(GeneName = Gene, SubjectID = Ind)


   # Gene-Ind level transformation
   cadd.collapse_grouped <- cadd.collapse_init %>% group_by(SubjectID, GeneName) %>%
     summarise(GC=max(GC),
            CpG=max(CpG),
            SIFTcat_deleterious=max(SIFTcat_deleterious),
            SIFTcat_tolerated=max(SIFTcat_tolerated),
            SIFTval=max(SIFTval),
            PolyPhenCat_benign=max(PolyPhenCat_benign),
            PolyPhenCat_possibly_damaging=max(PolyPhenCat_possibly_damaging),
            PolyPhenCat_unknown=max(PolyPhenCat_unknown),
            PolyPhenVal=max(PolyPhenVal),
            bStatistic=max(bStatistic),
            priPhCons=max(priPhCons),
            mamPhCons=max(mamPhCons),
            verPhCons=max(verPhCons),
            priPhyloP=max(priPhyloP),
            mamPhyloP=max(mamPhyloP),
            verPhyloP=max(verPhyloP),
            GerpN=max(GerpN),
            GerpS=max(GerpS),
            PHRED=max(PHRED),
            `EncodeH3K27ac-max`=max(`EncodeH3K27ac-max`), 
            `EncodeH3K4me1-max`=max(`EncodeH3K4me1-max`),
            `EncodeH3K4me3-max`=max(`EncodeH3K4me3-max`),
            `EncodeH3K27me3-max`=max(`EncodeH3K27me3-max`),
            `EncodeDNase-max`=max(`EncodeDNase-max`),
            `EncodeH3K36me3-max`=max(`EncodeH3K36me3-max`),
            `EncodeH3K9me3-max`=max(`EncodeH3K9me3-max`),
            `EncodeH3K9ac-max`=max(`EncodeH3K9ac-max`),
            `EncodeH3K4me2-max`=max(`EncodeH3K4me2-max`),
            `EncodeH3K79me2-max`=max(`EncodeH3K79me2-max`)) %>% 
  
   select(SubjectID, everything(.))
   #print(paste0("Collapsing Done: chr", chrom))

   # Impute missing values
   cadd.collapse_grouped = cadd.collapse_grouped %>%
     mutate(bStatistic = as.numeric(bStatistic)) %>% 
     mutate(GC=replace_na(GC, 0.418),
         CpG=replace_na(CpG, 0.024),
         SIFTcat_deleterious=replace_na(SIFTcat_deleterious,0),
         SIFTcat_tolerated=replace_na(SIFTcat_tolerated,0),
         SIFTval=replace_na(SIFTval, 0),
         PolyPhenCat_benign=replace_na(PolyPhenCat_benign, 0),
         PolyPhenCat_possibly_damaging=replace_na(PolyPhenCat_possibly_damaging, 0),
         PolyPhenCat_unknown=replace_na(PolyPhenCat_unknown, 0),
         PolyPhenVal=replace_na(PolyPhenVal,0),
         bStatistic=replace_na(bStatistic, 800.261),
         priPhCons=replace_na(priPhCons, 0.115),
         mamPhCons=replace_na(mamPhCons, 0.079),
         verPhCons=replace_na(verPhCons, 0.094),
         priPhyloP=replace_na(priPhyloP, -0.033),
         mamPhyloP=replace_na(mamPhyloP, -0.038),
         verPhyloP=replace_na(verPhyloP, 0.017),
         GerpN=replace_na(GerpN, 1.909),
         GerpS=replace_na(GerpS, -0.2),
         PHRED=replace_na(PHRED, 0),
        `EncodeH3K27ac-max`=replace_na(`EncodeH3K27ac-max`, 0), 
        `EncodeH3K4me1-max`=replace_na(`EncodeH3K4me1-max`, 0),
        `EncodeH3K4me3-max`=replace_na(`EncodeH3K4me3-max`, 0),
        `EncodeH3K27me3-max`=replace_na(`EncodeH3K27me3-max`, 0),
        `EncodeDNase-max`=replace_na(`EncodeDNase-max`, 0),
        `EncodeH3K36me3-max`=replace_na(`EncodeH3K36me3-max`, 0),
        `EncodeH3K9me3-max`=replace_na(`EncodeH3K9me3-max`, 0),
        `EncodeH3K9ac-max`=replace_na(`EncodeH3K9ac-max`, 0),
        `EncodeH3K4me2-max`=replace_na(`EncodeH3K4me2-max`, 0),
        `EncodeH3K79me2-max`=replace_na(`EncodeH3K79me2-max`, 0))

  #print(paste0("Gene Level Imputation Done: chr", chrom))

  # If it's the first chromosome, write with headers, otherwise append without headers
  if (chrom == 1) {
    fwrite(cadd.collapse_grouped, outfile1, quote = F, sep = '\t', row.names = F, col.names = TRUE)
  } else {
    fwrite(cadd.collapse_grouped, outfile1, quote = F, sep = '\t', row.names = F, col.names = FALSE, append = TRUE)
  }
  print(paste0("Writing of Gene level annotation done: chr", chrom))

  # copy gene-level annotation file to target annotation dir:
  #file.copy(outfile1, genelevel_dir)


  # Gene-Individual-Variant level transformation i.e, uncollapsed version
  cadd_filt <- cadd %>%
    group_by(Chrom, Pos, Ref, Alt) %>%
    arrange(desc(ConsScore)) %>%   # This will ensure max ConsScore is at the top for each group
    slice(1) %>%                    # Take the first row (which has the max value after arranging)
    ungroup()

  # Match rare variants with their annotations and rename columns
  cadd.collapse = filter(rv, Chrom == paste0("chr",chrom)) %>%
    left_join(., cadd_filt, by = c("Chrom" = "Chrom", "Start" = "Pos", "Ref"= "Ref", "Alt"="Alt")) %>%
    rename(GeneName = Gene, SubjectID = Ind)

  cadd.collapse$variantID <- paste(cadd.collapse$GeneName, paste(cadd.collapse$Chrom, cadd.collapse$End, cadd.collapse$Ref, cadd.collapse$Alt, sep = "_"), sep=":")

  # Derive the range of columns based on the names from cadd.collapse_grouped 
  # but apply this selection to the cadd.collapse dataframe
  # 1. Extract column names
  column_names <- names(cadd.collapse_grouped)

  # 2. Determine range of columns
  start_col <- which(column_names == "GC")
  end_col <- which(column_names == "EncodeH3K79me2-max")
  cols_to_select <- column_names[start_col:end_col]

  # 3. Select columns from cadd.collapse
  cadd.collapse_select <- cadd.collapse %>%
    select(SubjectID, GeneName, variantID, all_of(cols_to_select))

  print(paste0("Variant Level Concatenation Done: chr", chrom))

  # 4. Write selected Gene-Individual-Variant level data
  if (chrom == 1) {
    fwrite(cadd.collapse_select, outfile2, quote = F, sep = '\t', row.names = F, col.names = TRUE)
  } else {
    fwrite(cadd.collapse_select, outfile2, quote = F, sep = '\t', row.names = F, col.names = FALSE, append = TRUE)
  }

  print(paste0("Writing of Variant level annotation done: chr", chrom))

  # copy variant-level annotation file to target annotation dir:
  # file.copy(outfile2, variantlevel_dir)

}

# copy gene-level annotation file to target annotation dir:
file.copy(outfile1, genelevel_dir, overwrite = TRUE)

# copy variant-level annotation file to target annotation dir:
file.copy(outfile2, variantlevel_dir, overwrite = TRUE)

