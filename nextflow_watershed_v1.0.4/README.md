# Watershed Pipeline: Rare Genetic Variant Identification and Prioritization

---

## Table of Contents

1. [Watershed Pipeline Input Format Documentation](#watershed-pipeline-input-format-documentation)
2. [Other Required File Formats](#other-required-file-formats)
3. [Pipeline Overview & Strengths](#pipeline-overview--strengths)
4. [Installation & Environment Setup](#installation--environment-setup)
5. [Configuration & Profiles](#configuration--profiles)
6. [Running the Pipeline](#running-the-pipeline)
7. [Pipeline Steps & Outputs](#pipeline-steps--outputs)
8. [Enrichment and Outlier Methods](#enrichment-and-outlier-methods)
9. [Advanced Usage & Troubleshooting](#advanced-usage--troubleshooting)
10. [Watershed Model Details](#watershed-model-details)
11. [Example Workflows](#example-workflows)
12. [Contact](#contact)

---

## Watershed Pipeline Input Format Documentation

This document provides a detailed description of the input file formats required for the Nextflow watershed pipeline. Users must ensure that the input files conform to the specified formats to ensure proper execution of the pipeline. Below are the descriptions and formatting details for each of the five required input files:

### 1. TPM Input File (`tpm_infile`)

**Description:**  
This file contains Transcripts Per Million (TPM) data from RNAseq processing. The data input is similar to the input for typical eQTL analysis.

**Format:**  

- **File Type:** Tab-separated values (.tab)
- **Header:** The first row should contain sample names.
- **Columns:**
  - The first column should contain Ensembl gene IDs.
  - Subsequent columns should contain TPM values for each sample.

**Example:**

```
                    Sample1  Sample2  Sample3  Sample4
ENSG00000000003.15  0.123  0.456  0.789  0.012
ENSG00000000005.6   0.234  0.567  0.890  0.123
ENSG00000000419.14  0.345  0.678  0.901  0.234
ENSG00000000457.14  0.456  0.789  0.012  0.345
ENSG00000000460.17  0.567  0.890  0.123  0.456
ENSG00000000938.13  0.678  0.901  0.234  0.567
ENSG00000000971.16  0.789  0.012  0.345  0.678
ENSG00000001036.14  0.890  0.123  0.456  0.789
ENSG00000001084.13  0.901  0.234  0.567  0.890
```

### 2. Read Count Input File (`readcount_infile`)

**Description:**  
This file contains read count data from RNAseq processing, providing raw counts used in differential expression analysis.

**Format:**  

- **File Type:** Tab-separated values (.tab)
- **Header:** The first row should contain sample names.
- **Columns:**
  - The first column should contain Ensembl gene IDs.
  - Subsequent columns should contain read count values for each sample.

**Example:**

```
                    Sample1  Sample2  Sample3  Sample4
ENSG00000000003.15  10  20  30  40
ENSG00000000005.6   15  25  35  45
ENSG00000000419.14  20  30  40  50
ENSG00000000457.14  25  35  45  55
ENSG00000000460.17  30  40  50  60
ENSG00000000938.13  35  45  55  65
ENSG00000000971.16  40  50  60  70
ENSG00000001036.14  45  55  65  75
ENSG00000001084.13  50  60  70  80
```

### 3. Covariate Input File (`covariate_infile`)

**Description:**  
This file contains covariate data such as sex, principal components (PCs), and other potential confounders.

**Format:**  

- **File Type:** Tab-separated values (.tab/.tsv or .gz)
- **Header:** The first row should contain sample names.
- **Columns:**
  - The first column should contain covariate IDs.
  - Subsequent columns should contain covariate values for each sample.

**Example:**

```
id  Sample1  Sample2  Sample3  Sample4
sex  XY  XX  XY  XY
PC1  -0.001  -0.002  -0.003  -0.004
PC2  0.001  0.002  0.003  0.004
PC3  -0.005  -0.006  -0.007  -0.008
PC4  0.005  0.006  0.007  0.008
PC5  -0.009  -0.010  -0.011  -0.012
PEER1  -0.013  0.014  0.015  -0.016
PEER2  0.017  0.018  0.019  -0.020
PEER3  -0.021  -0.022  -0.023  0.024
PEER4  0.025  0.026  -0.027  0.028
PEER5  -0.029  -0.030  -0.031  -0.032
PEER6  0.033  -0.034  0.035  -0.036
PEER7  0.037  -0.038  0.039  0.040
PEER8  0.041  0.042  0.043  -0.044
```

### 4. Subject IDs File (`subjids_file`)

**Description:**  
This file contains the subject IDs along with their population ancestry.

**Format:**  

- **File Type:** Text file (.txt)
- **Columns:**
  - The first column should contain subject IDs.
  - The second column should contain a single population ancestry (e.g., EUR, AFR, EAS, SAS, GLOBAL).

**Example:**

```
Sample1  EUR
Sample2  EUR
Sample3  EUR
Sample4  EUR
Sample5  EUR
```

### 5. Rare Variant File (`rv_file`)

**Description:**  
This file contains rare variant data in VCF (Variant Call Format), including genotypic information for each sample.

**Format:**  

- **File Type:** VCF (.vcf.gz)
- **Header:** The VCF header lines start with `#` and include metadata. The last header line contains column names.
- **Columns:**
  - `#CHROM`: Chromosome number
  - `POS`: Position on the chromosome
  - `ID`: Variant identifier
  - `REF`: Reference allele
  - `ALT`: Alternate allele
  - `QUAL`: Quality score
  - `FILTER`: Filter status
  - `INFO`: Additional information
  - `FORMAT`: Format of the genotype fields
  - Sample columns containing genotypic data

**Example:**

```
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 Sample2 Sample3 Sample4
1       10000   rs0001  A       G       99.9    PASS    .       GT      0|0     0|1     1|1     0|0
1       20000   rs0002  C       T       99.8    PASS    .       GT      0|0     1|1     0|1     0|0
1       30000   rs0003  G       A       99.7    PASS    .       GT      1|1     0|0     0|1     1|1
1       40000   rs0004  T       C       99.6    PASS    .       GT      0|1     1|0     1|1     0|1
1       50000   rs0005  A       C       99.5    PASS    .       GT      0|0     0|0     1|0     1|1
```

---

## Other Required File Formats

- **Genomic Annotation Files:**  
  Used by the Watershed model, must be tab-separated with columns for SubjectID, GeneName, genomic annotations, outlier p-values, and N2pair.
- **YAML Environment Files:**  
  All environments are specified in `env_ymls/` (e.g., `watershed_env.yml`, `watershed_pyenv.yml`, etc).
- **Configuration Files:**  
  All Nextflow and resource configs are in `conf/` (e.g., `params.config`, `local.config`, etc).
- **Custom Scripts:**  
  All custom scripts and binaries are in `bin/` and are called by the pipeline as needed.

---

## Pipeline Overview & Strengths

The Watershed Nextflow pipeline is a **modular, scalable, and highly automated workflow** for rare variant identification and prioritization using RNA-seq and genomic data.  
**Key strengths and advantages:**

- **Automatic Sample Harmonization:**  
  The pipeline includes a robust "find common samples" feature. Users can provide sample lists from different input files (TPM, read counts, covariates, VCF, etc.), and the pipeline will automatically identify and harmonize the set of samples present in all required files. This ensures that downstream analyses are performed only on samples with complete data, reducing user error and manual preprocessing.

  **How it works:**  

  - You can provide input files with overlapping but not identical sample sets.
  - The pipeline will automatically detect the intersection of sample IDs across all files.
  - Only these common samples are used for all downstream steps.
  - This is handled internally, so users do not need to manually preprocess or match sample lists.

- **Flexible Input Handling:**  
  Accepts compressed and uncompressed files, supports various file extensions, and checks for format consistency.

- **Modular Design:**  
  Each major step (variant filtering, annotation, outlier detection, model training, etc.) is a separate module/subworkflow, making it easy to customize or extend.

- **Reproducibility:**  
  Uses Conda environments and supports Docker/Singularity for full reproducibility.

- **Scalability:**  
  Runs on local, SGE, SLURM, AWS Batch, and more. Resource profiles are easily customizable.

- **Comprehensive Reporting:**  
  Produces detailed output files, summary statistics, and publication-ready plots.

- **Advanced Variant Prioritization:**  
  Integrates the Watershed probabilistic model for rare variant effect prediction, and can support multi-phenotype outlier integration.

- **User-Friendly:**  
  Extensive documentation and example scripts.

---

## Installation & Environment Setup

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (>=20.10)
- [Conda](https://docs.conda.io/en/latest/) (with mamba recommended)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) (optional, for containers)
- R (>=3.5.1) and Python (>=3.6) for some scripts

### Environment Setup

All required environments are specified in `env_ymls/`. Example setup:

```bash
# Create the main environment
conda env create -f env_ymls/watershed_env.yml
conda activate watershed_env

# (Repeat for other YAMLs as needed)
```

---

## Configuration & Profiles

- **nextflow.config**: Main configuration file, includes all resource and parameter settings.
- **conf/**: Contains modular configs for different environments (local, SGE, SLURM, AWS, etc).
- **env_ymls/**: All Conda YAMLs for reproducible environments.

**Profiles:**  
Use `-profile local`, `-profile slurm`, `-profile sge`, `-profile conda`, etc., to select the appropriate compute environment.

---

## Running the Pipeline

### Basic Command

```bash
nextflow run ./main.nf -c nextflow.config \
  --tpm_infile "<path_to_tpm_infile>" \
  --readcount_infile "<path_to_readcount_infile>" \
  --covariate_infile "<path_to_covariate_infile>" \
  --subjids_file "<path_to_subjids_file>" \
  --rv_file "<path_to_rv_file>" \
  --analysis GLOBAL \
  --bin_dir "./bin" \
  --cache_dir "<cache_dir>" \
  --watershed_pyenv "env_ymls/watershed_pyenv.yml" \
  --watershed_renv "env_ymls/watershed_renv.yml" \
  --watershed_caddenv "env_ymls/watershed_caddenv.yml" \
  --watershed_env "env_ymls/watershed_env.yml" \
  --genotype_pcs 5 \
  --skip_cache \
  --cachedir_vep "<vep_cache_dir>" \
  -profile conda \
  -resume
```

### Comprehensive Example

Here's a complete example showing all available parameters with their typical values:

```bash
nextflow run ./main.nf -c nextflow.config \
  --cohort_name gtexEUR \
  --ancestry EUR \
  --outlier_method "pvalue" \
  --bin_dir "./bin" \
  --cache_dir "./cache" \
  --watershed_pyenv "env_ymls/watershed_pyenv.yml" \
  --watershed_renv "env_ymls/watershed_renv.yml" \
  --watershed_caddenv "env_ymls/watershed_caddenv.yml" \
  --watershed_env "env_ymls/watershed_env.yml" \
  --genotype_pcs 5 \
  --skip_cache \
  --cachedir_vep "$VEP_CACHEDIR" \
  --tpm_infile "<path_to_tpm_file>" \
  --readcount_infile "<path_to_readcount_file>" \
  --covariate_infile "<path_to_covariate_file>" \
  --subjids_file "<path_to_subjids_file>" \
  --rv_file "<path_to_rv_file>" \
  --tissue "Adipose" \
  --enrichment_method "both" \
  --aupr_plots "true" \
  -profile conda \
  -resume
```

**Replace** all `<...>` and `$...` with your actual file paths and environment variables.

### Parameters

- `--tpm_infile`: Path to TPM input file
- `--readcount_infile`: Path to read count input file
- `--covariate_infile`: Path to covariate input file
- `--subjids_file`: Path to subject IDs file
- `--rv_file`: Path to rare variant VCF file
- `--analysis`: Analysis type (e.g., GLOBAL)
- `--bin_dir`: Path to bin directory with scripts
- `--cache_dir`: Path to cache directory
- `--watershed_*env`: Paths to environment YAMLs
- `--genotype_pcs`: Number of genotype PCs to use
- `--skip_cache`: Skip cache step (optional)
- `--cachedir_vep`: Path to VEP cache directory
- `--cohort_name`: Name of the cohort being analyzed
- `--ancestry`: Population ancestry (e.g., EUR, AFR, EAS, SAS)
- `--outlier_method`: Method for outlier detection ("pvalue" or "zscore")
- `--enrichment_method`: Method for enrichment analysis ("pvalue", "zscore", or "both")
- `--aupr_plots`: Whether to generate AUPR plots ("true" or "false")
- `--tissue`: Name of the tissue being analyzed

---

## Pipeline Steps & Outputs

**Major Steps:**

1. **Sample Selection & Harmonization:**  
   - Automatically finds and uses only the set of samples present in all input files.
   - No need for manual sample matching.

2. **Rare Variant Processing:**  
   - Filters and annotates rare variants.

3. **Outlier Detection:**  
   - Identifies expression outliers.

4. **Annotation:**  
   - Annotates variants with VEP, CADD, and other tools.

5. **Annotation Collapsing:**  
   - Collapses and merges annotation data.

6. **N2 Pair Identification:**  
   - Finds N2 pairs for model evaluation.

7. **Watershed Model Run:**  
   - Runs the probabilistic model for variant prioritization.

8. **Visualization:**  
   - Generates enrichment and AUPR plots.

**Outputs:**

- Annotated variant files
- Outlier call files
- Merged annotation tables
- Watershed model results (RDS, posterior probabilities)
- Plots (enrichment, AUPR, etc.)

---

## Enrichment and Outlier Methods

The pipeline provides flexible options for both outlier detection and enrichment analysis, with detailed technical implementations:

### Outlier Detection Methods (`--outlier_method`)

- `"pvalue"`: Uses p-values to identify expression outliers
- `"zscore"`: Uses z-scores to identify expression outliers

**Description:**
The pipeline automatically selects the appropriate script based on the chosen outlier detection method.

### Enrichment Analysis Methods (`--enrichment_method`)

- `"pvalue"`: Uses p-value based enrichment analysis
- `"zscore"`: Uses z-score based enrichment analysis
- `"both"`: Runs both p-value and z-score based enrichment analyses

**Description:**
The pipeline uses different R scripts based on the method.

### AUPR Analysis (`--aupr_plots`)

- `"true"`: Generates Area Under Precision-Recall (AUPR) plots
- `"false"`: Skips AUPR plot generation

**Description:**
- The AUPR analysis provides comprehensive model evaluation through:
  - Precision-recall curves
  - Area under curve calculations
  - Model performance metrics

### Watershed Model Parameters

The model uses a sophisticated probabilistic framework that integrates multiple data types:

- **Genomic Annotations:** Incorporates various genomic features
- **RNA-seq Outlier Calls:** Integrates expression outlier data
- **Model Architecture:**
  - Conditional Random Field (CRF) model
  - Expectation-Maximization (EM) algorithm for parameter estimation
  - `--genotype_pcs`: Number of genotype principal components (default: 5)

### Model Training Parameters

The training process implements several advanced techniques:

- **Optimization:**
  - LBFGS optimization for logistic regression
  - Cross-validation for parameter selection
  - Dirichlet priors for parameter estimation

- **Features:**
  - Multiple dimensions for outlier detection
  - Flexible parameter initialization
  - Robust convergence criteria

### Output Generation

The pipeline produces comprehensive outputs:

- **Analysis Results:**
  - Enrichment plots (based on selected method)
  - AUPR plots (if enabled)
  - Model evaluation objects
  - Posterior probability files
  - Prediction objects

- **Model Artifacts:**
  - Trained model parameters
  - Performance metrics
  - Visualization files

### Example Configuration

```bash
nextflow run ./main.nf -c nextflow.config \
  --outlier_method "pvalue" \
  --enrichment_method "both" \
  --aupr_plots "true" \
  --genotype_pcs 5 \
  # ... other parameters ...
```

**Note:** The `"both"` option for enrichment method is particularly useful for comparing results from different approaches and validating findings across methods.

---

## Advanced Usage & Troubleshooting

- **Resuming:** Use `-resume` to continue a previous run.
- **Custom Modules:** Add or modify subworkflows in `modules/` and `subworkflows/`.
- **Resource Tuning:** Edit `conf/profile_resources.config` for custom resource allocation.
- **Troubleshooting:** Check `.nextflow.log`.
- **Flexible Input:** Accepts compressed/uncompressed files, various extensions, and checks for format consistency.

---

## Watershed Model Details

Watershed is an unsupervised probabilistic framework that integrates genomic annotations and RNA-seq outlier calls to estimate the probability that a rare variant has a functional effect.

**Input Format:**

- Each line: gene-individual pair
- Columns:
  1. SubjectID
  2. GeneName
  3. Genomic annotation columns (user-defined)
  4. Outlier p-value columns (user-defined, can be signed for direction)
  5. N2pair (unique integer for N2 pairs, "NA" otherwise)

**Example:**

```
SubjectID GeneName annot1 annot2 ... outlier1 outlier2 ... N2pair
```

**Key Scripts:**

- `evaluate_watershed.R`: Trains and evaluates the model.
- `predict_watershed.R`: Predicts posterior probabilities on new data.

**Dependencies:**

- R 3.5.1, Rcpp, lbfgs, PRROC, optparse

See `src/Watershed/README.md` for full details.

---

## Example Workflows

### Minimal Example

```bash
nextflow run ./main.nf -c nextflow.config \
  --tpm_infile data/tpm.tab \
  --readcount_infile data/readcount.tab \
  --covariate_infile data/covariates.tab \
  --subjids_file data/subjects.txt \
  --rv_file data/rare_variants.vcf.gz \
  -profile local \
  -resume
```

### Using SLURM

```bash
nextflow run ./main.nf -c nextflow.config \
  ... (other params) ... \
  -profile slurm \
  -resume
```

### Customizing Resources

Edit `conf/profile_resources.config` and select with `-profile` as needed.

---

## Contact

**Pipeline Author:**  

- Surya B. Chhetri

**Contact:**  

- chhetribsurya@gmail.com

---

**For more details, see the documentation in each subfolder.**

---