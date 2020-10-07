#!/usr/bin/env Rscript

# Script to take base GWAS catalogue summary statistics and format them so that they match SAIGE statistics and can be used for input to PRSice input:

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

args= commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("this script requires the following input: \n A file containing GWAS catalogue summary statistics (.tsv)")
}

gwas_catalogue_file <- as_tibble(fread(args[1]))

# Transform and clean up data

base <- gwas_catalogue_file %>% select(starts_with("hm_"), "p-value")

base <- filter(base, !(is.na(hm_beta) == TRUE & is.na(hm_odds_ratio) == TRUE))
