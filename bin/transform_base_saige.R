#!/usr/bin/env Rscript

# Script to take base saige output (obtained from lifebit-ai/gel-gwas) and format it for input to PRSice input:

suppressPackageStartupMessages({
library(tidyverse)
library(data.table)
})

args= commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("this script requires the following input: \n a SAIGE file (.csv) obtained from lifebit-ai/gel-gwas")
}

saige_file <- as_tibble(fread(args[1]))

base <- saige_file
# TODO: Add some sanity checks:
# - Any SNPs with no betas/ORs?

write.table(base, "base.data", quote = F, row.names =F, sep = " ")


