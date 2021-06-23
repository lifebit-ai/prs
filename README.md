# Post-GWAS polygenic risk scores (PRS) pipeline

## 1 - Description

This pipeline calculates polygenic risk scores (PRS) using `PRSice`. It requires the following inputs:

- **A target cohort**: genotype data of the individuals for which you wish to obtain PRS scores. 

- **A target phenotype**: a file specifying sample IDs, phenotypes and covariates.

- **A base cohort**: GWAS summary statistics which will be used to calculate the PRS. These can be either:
  
  - SAIGE formatted summary statistics (i.e. an output of the `lifebit-ai/biobank-gwas` pipeline)
  
  - GWAS catalogue summary statistics (still in development)

For more usage details, including a description of all available parameters, see `docs/usage.md`

## 2 - Important assumptions

In its current implementation, the pipeline assumes it is being run post-GWAS, more specifically, after having run the `lifebit-ai/biobank-gwas` **via the CloudOS Cohort Browser**. It therefore assumes the following about the target dataset and its correspond phenotype file:

- This PLINK files are the output of the `lifebit-ai/biobank-gwas` pipeline. These PLINK files have therefore already undergone standard QC.

- In accordance with the previous point, the phenotype file must correspond to the phenotype file used by `lifebit-ai/biobank-gwas` (run via the CloudOS Cohort Browser). See example in `testdata/cohort_parsed_file.phe`.

- The target data consists of PLINK files, split by chromosome. See example in `testdata/plink`.

##  3 - Example usage (using SAIGE summary statistics as base cohort):

```
nextflow run main.nf \
--saige_base testdata/saige_results_covid_1.csv \
--target_pheno testdata/cohort_parsed_file.phe \
--target_plink_dir testdata/plink
```

For more usage details, including a description of all available parameters, see `docs/usage.md`.

