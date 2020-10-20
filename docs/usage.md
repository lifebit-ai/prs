# Usage documentation

# 1 - Introduction

This document provides a detailed guide on how to run the `lifebit-ai/prs` pipeline. This pipeline uses to the `PRSice` package to calculate polygenic risk scores (PRS).

# Usage

In order to calculate a PRS, one must provide:

## (1) a target cohort

This corresponds to the genotype data of the individuals for which you wish to obtain PRS scores. In its current implementation, the pipeline assumes it is being run post-GWAS, more specifically, after having run the `lifebit-ai/gel-gwas` **via the CloudOS Cohort Browser**. It therefore assumes the following about the target dataset and its correspond phenotype file:

- This PLINK files are the output of the `lifebit-ai/gel-gwas` pipeline. These PLINK files have therefore already undergone standard QC.

- The target data consists of PLINK files, split by chromosome. See example in `testdata/plink`

## (2) a target cohort phenotype file

This corresponds to a a file specifying sample IDs, phenotypes and covariates. See example in `testdata/cohort_parsed_file.phe`. In accordance with the previous section, the phenotype file must correspond to the phenotype file used by `lifebit-ai/gel-gwas` (run via the CloudOS Cohort Browser).

## (3) base cohort

This file contains GWAS summary statistics which will be used to calculate the PRS. These can be either:

   - **SAIGE summary statistics**
 
     Obtained from a previously run GWAS. This file should have a `.csv` format. This can be an output of the `lifebit-ai/gel-gwas` pipeline.

   - **GWAS catalogue summary statistics**

     Obtained from the GWAS catalogue. This file should have a `.h.tsg.gz` format.

NB: Only 1 base can be used by the pipeline.

# 2- Basic example

- When using SAIGE summary statistics as your base cohort:

```
nextflow run main.nf \
--saige_base testdata/saige_results_covid_1.csv \
--target_pheno testdata/cohort_parsed_file.phe \
--target_plink_dir testdata/plink
```

- When using a GWAS catalogue summary statistics as your base cohort:

```
nextflow run main.nf \ 
--gwas_cat_study_id GCST004420 \
--target_pheno testdata/cohort_parsed_file.phe \
--target_plink_dir testdata/simulated-data/updated-simulated-plink/updated-split
```

Where `GCST004420` represents the study accession ID of a GWAS reported in the GWAS catalogue. See: https://www.ebi.ac.uk/gwas/studies/GCST004420
Note that special simulated PLINK files were made for the purpose of this test command. See more about this in `testdata/simulated-data/updated-simulated-plink/updated-split`

# 3- Essential parameters

**--saige_base**: the SAIGE summary statistics (`.csv`) file that will be used as a base cohort.

**--gwas_cat_study_id**: the study accession ID of a GWAS reported in the GWAS catalogue. See: https://www.ebi.ac.uk/gwas/studies/GCST004420

**--target_pheno**: a phenotype file contain sample IDs, phenotype IDs and covariates.

**--target_plink_dir**: path to folder containing PLINK files of the target cohort.

# 4- Other parameters

**--binary_trait**: whether the phenotype is binary (`true`) or continous (`false`). Default: `true`. 


