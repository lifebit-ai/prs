#!/usr/bin/env Rscript

# Script to take base data (SAIGE or GWAS catalogue), new coordinates and update the build of the base data



###################
# Import packages #
###################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(optparse)
})



#####################
# Parsing arguments #
#####################

option_list = list(
  make_option(c("--input_base"), action="store", type='character', help="String containing input base.data file"),
  make_option(c("--input_coordinates"), action="store", type='character', help="String containing input coordinates file (i.e. the new build)"))

args = parse_args(OptionParser(option_list = option_list))

# Args to variable
input_base              = args$input_base
input_coordinates       = args$input_coordinates



######################################
# Importing data and transforming it #
######################################

base <- as_tibble(fread(input_base))
coord <- as_tibble(fread(input_coordinates))

coord$Allele1 <- str_sub(coord$genotype, 1,1)
coord$Allele2 <- str_sub(coord$genotype, 2,2)

base$CHR <- coord$chromosome
base$POS <- coord$position
base$Allele1 <- coord$Allele1
base$Allele2 <- coord$Allele2


write.table(base, "matching.base.data", quote = F, row.names =F, sep = " ")



