#!/usr/bin/env Rscript

# Script to take base saige output (obtained from lifebit-ai/biobank-gwas) and format it for input to PRSice input:



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

#### Remove SNPs with no beta - these cannot be used by PRSice ####

base <- filter(saige_file, !(is.na(BETA) == TRUE))

#### For SNPs with no p-value, replace the NA by a 1 ####

base <- base %>% mutate(p.value = if_else(is.na(p.value), 1, p.value))

#### Remove duplicate SNPs - these cannot be used by PRSice (an error will be thrown) ####

base <- distinct(base, SNPID, .keep_all = TRUE)

#### Save data ####

write.table(base, "base.data", quote = F, row.names =F, sep = " ")


