#!/usr/bin/Rscript --vanilla
library(readr)
library(Matrix)

tsv_Files <- commandArgs(trailingOnly = TRUE)

rds_Files <- gsub(pattern = ".tsv.gz$", replacement = ".rds", x = tsv_Files)

Import_umi_tools_Gene_x_CelL_Matrix <- function(path_to_tsv){
  ###Import tsv.gz object
  dense_matrix <- read_tsv(path_to_tsv)
  ###Isolate gene names
  genes <- dense_matrix$gene
  ###Remove gene names from data
  dense_matrix <- dense_matrix[,2:ncol(dense_matrix)]
  ###Isolate cell names
  cells <- colnames(dense_matrix)
  ###Convert from tibble to dense matrix
  dense_matrix <- as.matrix(dense_matrix)
  ###Convert from dense matrix to sparse matrix
  sparse_matrix <- Matrix::Matrix(data = dense_matrix, sparse = TRUE,
                                  dimnames = list(genes, cells))
  ###Return object
  return(sparse_matrix)
}

matrices <- Import_umi_tools_Gene_x_CelL_Matrix(tsv_Files)

saveRDS(matrices, file = rds_Files)
