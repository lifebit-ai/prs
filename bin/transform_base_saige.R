#!/usr/bin/env Rscript

# Script to take base saige output (obtained from lifebit-ai/gel-gwas) and format it for input to PRSice input:



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

option_list <- list(make_option(c("--input_saige"), action="store", type='character',help="String containing input SAIGE file (base cohort)"))
args = parse_args(OptionParser(option_list = option_list))

# Arg to variable
input_saige = args$input_saige



######################################
# Importing data and transforming it #
######################################

saige_file <- as_tibble(fread(input_saige))

base <- saige_file
# TODO: Add some sanity checks:
# - Any SNPs with no betas/ORs?

write.table(base, "base.data", quote = F, row.names =F, sep = " ")


