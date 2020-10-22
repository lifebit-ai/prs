#!/usr/bin/env Rscript

# Script to plot the relationship between the "best-fit" PRS and thhe covariates used for the GWAS and PRS



######################
# Importing packages #
######################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(snakecase)
  library(optparse)
})



#####################
# Parsing arguments #
#####################

option_list = list(
  make_option(c("--input_cov"), action="store", type='character', help="String containing input cov file"),
  make_option(c("--input_prs"), action="store", type='character', help="String containing input prs results file (PRSice.best)"),
  make_option(c("--input_metadata"), action="store", type='character', help="String containing input phenotype metadata file"))

args = parse_args(OptionParser(option_list = option_list))

# Args to variable
input_cov       = args$input_cov
input_prs       = args$input_prs
input_metadata  = args$input_metadata



######################################
# Importing data and transforming it #
######################################

cov <- as_tibble(fread(input_cov))
prs <- as_tibble(fread(input_prs))
pheno_metadata <- as_tibble(fread(input_metadata))

#### Used pheno metadata to keep only relevant covariates for plotting: continuous and/or numeric variables ####

colnames(pheno_metadata) <- snakecase::to_snake_case(colnames(pheno_metadata))
pheno_metadata$field_name <- snakecase::to_snake_case(pheno_metadata$field_name)

colnames(cov) <- snakecase::to_snake_case(colnames(cov))
colnames(prs) <- snakecase::to_snake_case(colnames(prs))

# Obtain covariates of interest and remove any year-base ones: year of birth for example

cov_to_keep <- pheno_metadata %>% filter(field_id_type == "Continuous" | field_id_type == "Integer") %>% select(field_name)
year_based_cov <- str_which(cov_to_keep$field_name, "year")
cov_to_keep <- cov_to_keep[-year_based_cov,]
cov_to_keep <- as.character(cov_to_keep$field_name)

plot_data <- inner_join(cov, prs, by = c( "iid" = "iid", "fid" = "fid"))
plot_data <- suppressWarnings(plot_data %>% select(iid, prs, one_of(cov_to_keep)))

#### Plot ####

plot_data <- as.data.frame(plot_data)

suppressWarnings(
for (i in colnames(plot_data[, -which(names(plot_data) %in% c("iid", "prs")), drop = F])) {
  
  cov_name <- i 
  
  plot <- ggplot(plot_data, aes(x=prs, y=plot_data[[cov_name]]))+
    geom_point()+
    theme_classic()+
    labs(x="Polygenic Risk Score", y=cov_name)
  
  plot.saved<-ggsave(plot, file= paste("prs_vs_",cov_name,".png", sep = ""))
  
}
)


