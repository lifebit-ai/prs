# Making simulated PLINK data

## 1 - Aim

The purpose of the following is to simulate PLINK genotype files such that:
- The PLINK files have enough SNPS and samples for testing and match the numbers I have in my test phenotype file (`cohort_parsed_file.phe`)
- The PLINK files are split per chromosome
- The PLINK files contain SNPs that will match those found in a GWAS catalogue summary statistics file: this way the pipeline can be appropriately tested.

For the purpose of this testing, the following GWAS catalogue study was chosen (study GCST004420):
```
$ cd testdata
$ mkdir gwas-catalogue-data && cd gwas-catalogue-data
$ wget ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/Ahola-OlliAV_27989323_GCST004420/harmonised/27989323-GCST004420-EFO_0008082.h.tsv.gz
$ gunzip -d 27989323-GCST004420-EFO_0008082.h.tsv.gz
```

Due to the fact this file is very large (1.2GB), it is not commited to this repository.

## 2 - Software

In order to simulate PLINK files, one can use PLINK1:
https://zzz.bwh.harvard.edu/plink/simulate.shtml

Important considerations:
- All SNPs simulated are unlinked and in linkage equilibrium.
- This tool only generates individuals drawn from a homogeneous population (but see link for tip on how to get around that)
- You can also use PLINK to simulate summary statistics directly

## 3 - Commands used

To produce additional test for this pipeline, the following commands were used:

```
#### 1 -  Obtain PLINK1 docker image ####

$ docker pull lifebitai/plink1:latest

#### 2 -  Make a wgas.sim file and simulate data ####

# Here, we want to simulate 2200 SNPs (100 per chromosome) for 112 individuals (7 cases, 105 controls)
# These numbers should make it easier for me to split the data and use the current pheno file I have used for testing the pipeline.

# Use the Docker image to easily run PLINK
$ pwd
# [...]testdata/simulated-data
$ docker run --rm -it -v "$PWD":"$PWD" -w "$PWD" -e HOME=$PWD” --user "$(id -u):$(id -g)" lifebitai/plink1:latest

$ plink --simulate original-simulated-plink/wgas.sim --simulate-ncases 7 --simulate-ncontrols 105 --out original-simulated-plink/simdata
$ plink --bfile original-simulated-plink/simdata --recode --out original-simulated-plink/simdata

# IMPORTANTLY: 
# - (1) The pheno column present in the final version of the .fam files (split per chromosome will have -9 regardless of being cases and controls)
#   - Indeed, our external pheno file will take care of handling phenos (just like it has when using SAIGE statistics)
#   - The reason we still used --simulate-ncontrols and --simulate-ncases is because it seems to be the only way of controlling the final number of participants.
# - (2) I get an error about some options not being compatible with --simulate when I used --recode in the same command so I split the command in 2 steps.

#### 3 - Use R script to update simulated data to obtain matching SNPS with GWAS catalogue data ####

# Run update_simdata.R
# NB: the Docker image did not have R so I ran the interactively via Rstudio (possible due to the way the image is mounted)

#### 4 - Use PLINK to update the simdata a final time ####

# First, ensure seperations are correct
$ cat updated-simulated-plink/update_allele.txt | awk '{print $1"\t"$2 " " $3"\t"$4 " " $5}' > updated-simulated-plink/update_allele-modified.txt

$ plink --file updated-simulated-plink/simdata --update-alleles updated-simulated-plink/update_allele-modified.txt --allow-no-sex --out updated-simulated-plink/simdata-final --make-bed

#### 5 - Use PLINK to split per chromosome and output bed files ####

for i in {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}
do
    plink \
    --bfile updated-simulated-plink/simdata-final \
    --chr $i \
    --allow-no-sex \
    --out updated-simulated-plink/updated-split/sampleA_chr${i}_filtered \
    --make-bed
done
```

# - 4 Making a subset of GWAS catalogue data in order to make a CI testing profile

The aim here was to make a subset of the a GWAS catalogue that we can use for making a CI testing profile that tests the ability of the pipeline to handle GWAS catalogue data.

```
$ cd testdata/simulated-data/gwas-catologue-data
$ head -n 350000 27989323-GCST004420-EFO_0008082.h.tsv > subset_GCST004420.h.tsv
```

Next steps:
- This file was committed to GitHub. 
- Then it a download link to its GitHub location was added to `assets/ftp_locations_harmonized.csv` in order to perform a download.

Benefits:
1) The pipeline can then be run as a test locally, without needing to use 16CPUs for the `transform_gwas_catalogue_base` step

```
nextflow run main.nf \
--gwas_cat_study_id GCST004420-ci \
--target_pheno testdata/cohort_parsed_file.phe \
--target_plink_dir testdata/simulated-data/updated-simulated-plink/updated-split
```

2) This subset GWAS catalogue data is used for CI testing (indeed, one cannot use a real size dataset for CI testing as these are too large).


