#!/usr/bin/env Rscript

# Script to take target pheno (obtained from lifebit-ai/gel-gwas) and format it for input to PRSice input:
# - .pheno file
# - .cov file 

suppressPackageStartupMessages({
library(tidyverse)
library(data.table)
})

args= commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("this script requires the following input: \n a pheno file (.phe) obtained from lifebit-ai/gel-gwas")
}

phe_file <- as_tibble(fread(args[1]))

pheno <- phe_file  %>% select("FID", "IID", "PHE")
cov <- phe_file  %>% select(-"PHE")

write.table(pheno, "target.pheno", quote = F, row.names =F, sep = " ")
write.table(cov, "target.cov", quote = F, row.names =F, sep = " ")


