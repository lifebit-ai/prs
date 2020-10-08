#!/usr/bin/env Rscript

# Script to plot to make a PRS density plot of cases vs controls

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(snakecase)
  library(sm)
})

args= commandArgs(trailingOnly=TRUE)

#### Function ####

plot_density <- function(prs.data) {
  
  # (1) Create vector of prs scores for cases/controls
  cases = prs.data[(prs.data$phe == "1"),]$prs
  controls = prs.data[(prs.data$phe == "0"),]$prs
  
  # (2) Create a vector of groups labels - one observation = group sample belongs to
  group.index <- rep(1:2, c(length(cases), length(controls)))
  
  # (3) Use sm.density.compare() function to calculate density for all groups and plot them
  den = sm.density.compare(c(cases,controls), group = group.index, model = "none", xlab = "PRS scores",lty = c(1,1,1,1))
  
  # (4) Add legend to plot
  legend("topright", legend = c("cases","controls"), fill = 2+(0:4),cex = 0.9)
  
}

#### Import data ####

pheno <- as_tibble(fread(args[1]))
prs <- as_tibble(fread(args[2]))

#### Transform data ####

colnames(pheno) <- snakecase::to_snake_case(colnames(pheno))
colnames(prs) <- snakecase::to_snake_case(colnames(prs))
plot_data <- inner_join(pheno, prs, by = c( "iid" = "iid", "fid" = "fid"))

#### Plot density ####

png("prs_vs_density.png")
plot_density(plot_data)
dev.off()


