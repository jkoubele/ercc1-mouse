library(Seurat)
library(tidyverse)
library(HGNChelper)
library(openxlsx)

source("/mnt/ercc1/ercc1-mouse/scripts/utils.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R") 
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")


path_prefix <- "/mnt/ercc1"

input_folder <- file.path(path_prefix, "clustering_output")
output_folder <- file.path(path_prefix, "cell_type_annotation")

seurat_object <- readRDS(file.path(input_folder, "seurat_object_after_clustering.rds"))


# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
tissue <- "Kidney"

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

system(paste("chmod -R 777", output_folder))
