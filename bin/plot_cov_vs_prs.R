#!/usr/bin/env Rscript

# Script to plot the relationship between the "best-fit" PRS and thhe covariates used for the GWAS and PRS

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(snakecase)
})

args= commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("this script requires the following inputs: \n - the pheno file used for PRSice, \n - the covariate file used for PRSice \n - the PRSice.best file produced by PRSice")
}

#### Import data ####

pheno <- as_tibble(fread(args[1]))
cov <- as_tibble(fread(args[2]))
prs <- as_tibble(fread(args[3]))

pheno_metadata <- as_tibble(fread("assets/Metadata phenotypes - Mapping file.csv"))

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
for (i in colnames(plot_data[-c(1,2)])) {
  
  cov_name <- i 
  
  plot <- ggplot(plot_data, aes(x=prs, y=plot_data[[cov_name]]))+
    geom_point()+
    theme_classic()+
    labs(x="Polygenic Risk Score", y=cov_name)
  
  plot.saved<-ggsave(plot, file= paste(cov_name,"_vs_prs.png", sep = ""))
  
}
)


