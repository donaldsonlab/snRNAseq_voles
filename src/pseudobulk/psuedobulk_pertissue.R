#! /usr/bin/env Rscript

#-------------------------------------------------------------------------------
# Script Name: psuedobulk_pertissue.R
# Purpose: Create output files that are per celltype counts per animal.
# Author: Mary Allen
# Date: July 28, 2025
#-------------------------------------------------------------------------------


library(Seurat)
library(data.table)
library(DESeq2)
library(tidyverse)


#Set paths
outdir= "/Users/allenma/snRNAseq_voles_files/output/" 
aggr_counts <- read.csv("/Users/allenma/snRNAseq_voles_files/output/psuedobulk.csv")


