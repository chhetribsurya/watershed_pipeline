#!/usr/bin/env bash

# load R
ml r/4.2.0

# RUN Collapse annotation scripts
Rscript step12.1_collapse_vep_loftee_Final.R
Rscript step12.2_collapse_cadd.R
Rscript step12.3_collapse_ucsc.R
Rscript step12.4_collapse_gencodeDistFromGene.R
Rscript step12.5_collapse_kgpexFreq.R


