#!/usr/bin/env Rscript

# Script to take base GWAS catalogue summary statistics and format them so that they match SAIGE statistics and can be used for input to PRSice input:



###################
# Import packages #
###################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args= commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("this script requires the following input: \n A file containing GWAS catalogue summary statistics (.tsv)")
}



###############
# Import data #
###############

gwas_catalogue_file <- as_tibble(fread(args[1]))



###############################
# Transform and clean up data #
###############################

#### Keep harmonized data only ####

base <- gwas_catalogue_file %>% select(starts_with("hm_"), starts_with("p_"))

#### Remove SNPs with no beta or OR - these cannot be used by PRSice ####

base <- filter(base, !(is.na(hm_beta) == TRUE & is.na(hm_odds_ratio) == TRUE))

#### For SNPs with at least a beta or OR, alternatively use the beta or OR to calculate the other ####

base$hm_beta <- as.numeric(base$hm_beta)
base$hm_odds_ratio <- as.numeric(base$hm_odds_ratio)

base <- base %>%
  mutate(hm_beta = if_else(is.na(hm_beta), log(hm_odds_ratio), hm_beta), 
         hm_odds_ratio = if_else(is.na(hm_odds_ratio), exp(hm_beta), hm_odds_ratio))

#### Remove SNPs with no p-value - these cannot be used by PRSice (which needs to select SNPs based on P-value thresholds) ####

base <- filter(base, !(is.na(p_value) == TRUE))
  
#### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####

base <- distinct(base, hm_rsid, .keep_all = TRUE)
base <- distinct(base, hm_variant_id, .keep_all = TRUE)

#### Change column names to match SAIGE output ####

colnames(base) <- c("hm_variant_id", "SNPID", "CHR", "POS", "Allele2", "Allele1", "BETA", "OR", "hm_ci_lower","hm_ci_upper","hm_effect_allele_frequency","hm_code","p.value")

#### Save data ####

write.table(base, "base.data", quote = F, row.names =F, sep = " ")


