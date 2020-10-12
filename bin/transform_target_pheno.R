#!/usr/bin/env Rscript

# Script to take target pheno (obtained from lifebit-ai/gel-gwas) and format it for input to PRSice input:
# - .pheno file
# - .cov file 



######################
# Importing packages #
######################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})



#####################
# Parsing arguments #
#####################

option_list <- list(make_option(c("--input_pheno"), action="store", type='character',help="String containing input pheno file"))
args = parse_args(OptionParser(option_list = option_list))

# Arg to variable
input_pheno = args$input_pheno



######################################
# Importing data and transforming it #
######################################

phe_file <- as_tibble(fread(input_pheno))

pheno <- phe_file  %>% select("FID", "IID", "PHE")
cov <- phe_file  %>% select(-"PHE")

# Remove any covariate which has at least 1 NA: indeed PRSice does not seem to like having covariates with NAs and will throw an error. 
# For now, remove the relevant covariates but a card has been added to repo to better deal with this in the future.
cov <- cov %>% select(where(~!any(is.na(.))))

write.table(pheno, "target.pheno", quote = F, row.names =F, sep = " ")
write.table(cov, "target.cov", quote = F, row.names =F, sep = " ")


