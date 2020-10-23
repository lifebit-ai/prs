#!/usr/bin/env Rscript

# Script to update simulated PLINK file so that the SNPs etc match thoses from a GWAS catalogue:
# - Takes simulated PLINK data
# - Takes in GWAS catalogue data
# - Takes in pheno file



###################
# Import packages #
###################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})



#####################################################
# Importing GWAS catalogue data and transforming it #
#####################################################

gwas_catalogue_file <- as_tibble(fread("gwas-catalogue-data/27989323-GCST004420-EFO_0008082.h.tsv"))

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

#### Obtain exactly 2200 SNPs - 100 per chromosome ####

# First remove rows with X chromosome (won't use them for my simulated data)
base <- base %>% filter(!(CHR == "X"))

# Group by chromosome and select 100 SNPs for each
base_2200 <- base %>%
  group_by(CHR) %>%
  slice_sample(n = 100)



###########################################
# Importing simulated data and pheno file #
###########################################

ped <- as_tibble(fread("original-simulated-plink/simdata.ped"))
map <- as_tibble(fread("original-simulated-plink/simdata.map"))

pheno <- as_tibble(fread("../cohort_parsed_file.phe"))



###################################################
# Update PLINK files based on GWAS catalogue data #
###################################################

#### First, order base_2200 by chromosome and position ####

base_2200 <- base_2200[order(base_2200$CHR, base_2200$POS),] 

#### Update map ####

updated_map <- map
updated_map$V1 <- base_2200$CHR
updated_map$V2 <- base_2200$SNPID
updated_map$V4 <- base_2200$POS
#NB: no need to update V3 (can all be 0)

write.table(updated_map,"updated-simulated-plink/simdata.map", col.names = F, row.names = F, quote = F, sep = "\t")

#### Updating ped ##### 

updated_ped <- ped
updated_ped$V1 <- pheno$FID
updated_ped$V2 <- pheno$IID
updated_ped$V3 <- pheno$PAT
updated_ped$V4 <- pheno$MAT
updated_ped$V5 <- pheno$SEX
updated_ped$V6 <- -9

write.table(updated_ped,"updated-simulated-plink/simdata.ped", col.names = F, row.names = F, quote = F, sep = " ")

#### Update allele encoding: current only using d and D ####

# To do this, we will use PLINK's update-allele-information.
# So here, we just just need to make a file we will use with PLINK.

update_allele <- base_2200 %>% ungroup() %>% select("SNPID", "Allele1", "Allele2")
update_allele <- update_allele %>% mutate(old_allele1 = "d", old_allele2 = "D") %>% relocate("SNPID","old_allele1","old_allele2","Allele1","Allele2" )

write.table(update_allele,"updated-simulated-plink/update_allele.txt", col.names = F, row.names = F, quote = F, sep = " ")
