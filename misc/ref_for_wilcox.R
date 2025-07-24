rm(list=ls())

## Setup libraries -----------------------------------------------------------
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)

## function -----------------------------------------------------------
# Add decimal suffix to duplicate IDs to ensure uniqueness
add_decimal_numbers <- function(target_values, value_list) {
  target_indices <- match(value_list, target_values)
  modified_list <- ifelse(
    !is.na(target_indices),
    paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
    value_list
  )
  return(modified_list)
}

## Data import ---------------------------------------------------------------
genedata_org_hg18 <- fread("allcount_matrix_2023-04-11.txt", data.table=FALSE)
rownames(genedata_org_hg18) <- genedata_org_hg18$Geneid
genedata_org_hg18$Geneid <- NULL

genedata_org_hg19 <- fread("allcount_matrix_2023-04-12.txt", data.table=FALSE)
rownames(genedata_org_hg19) <- genedata_org_hg19$Geneid
genedata_org_hg19$Geneid <- NULL

genedata_org_v30  <- fread("allcount_matrix_2023-04-11.txt", data.table=FALSE)
rownames(genedata_org_v30)  <- genedata_org_v30$Geneid
genedata_org_v30$Geneid  <- NULL

genedata_org_v43  <- fread("allcount_matrix_2023-04-11.txt", data.table=FALSE)
rownames(genedata_org_v43)  <- genedata_org_v43$Geneid
genedata_org_v43$Geneid  <- NULL

genedata_org_T2T  <- fread("allcount_matrix_2023-07-03.txt", data.table=FALSE)

map               <- fread("map_genename_T2T_ensembl_from_gtf.txt")
assembly_map      <- readRDS("rbind_map_between_assemblies_10.04.23.RDS")

## Prepare assembly grid -----------------------------------------------------
v1 <- c("genedata_org_hg18","genedata_org_hg19",
        "genedata_org_v30","genedata_org_v43","genedata_org_T2T")
v2 <- c("hg18","hg19","v30","v43","T2T")
grid <- as.data.frame(cbind(v1, v2))

## Process each assembly -----------------------------------------------------
# Filter expressed genes, rename IDs, and convert to matrix form
for (x in seq_len(nrow(grid))) {
  print(grid[x,2])
  if (grid[x,2] != "T2T") {
    og_ref <- get(grid[x,1])
    ref    <- og_ref[,6:644]
    ref_isexpr <- rowSums(cpm(ref) >= 1) >= 0.1 * ncol(ref)
    ref    <- ref[ref_isexpr, ]
    ref$Geneid_noVersion <- unlist(lapply(
      strsplit(rownames(ref), ".", fixed=TRUE),
      function(x) {
        if (!any(grepl("_", x))) x[1] else paste(x[1], paste(strsplit(x[2], "_")[[1]][2:3], collapse = "_"), sep = "_")
      }
    ))
    rownames(ref) <- ref$Geneid_noVersion
    ref$Geneid_noVersion <- NULL
    ref <- as.matrix(ref)
    assign(paste0("matrix", grid[x,2]), ref)
  } else {
    og_ref <- get(grid[x,1])
    T2T    <- merge(og_ref, map, by.x = "Geneid", by.y = "gene_id")
    T2T$Geneid_noVersion <- unlist(lapply(
      strsplit(T2T$projection_parent_gene, ".", fixed=TRUE),
      function(x) {
        if (!any(grepl("_", x))) x[1] else paste(x[1], paste(strsplit(x[2], "_")[[1]][2:3], collapse = "_"), sep = "_")
      }
    ))
    dupl <- T2T[duplicated(T2T$projection_parent_gene), ]
    uni_ids <- unique(dupl$projection_parent_gene)
    uni_ids <- data.frame(ids = uni_ids,
                          Geneid_noVersion = unlist(lapply(strsplit(uni_ids, ".", fixed=TRUE), `[`, 1)))
    modified_values <- add_decimal_numbers(uni_ids$Geneid_noVersion, T2T$Geneid_noVersion)
    rownames(T2T) <- modified_values
    T2T$projection_parent_gene <- NULL
    T2T$Geneid_noVersion <- NULL
    T2T_new <- T2T[,7:645]
    T2T_new_isexpr <- rowSums(cpm(T2T_new) >= 1) >= 0.1 * ncol(T2T_new)
    T2T_new <- T2T_new[T2T_new_isexpr, ]
    matrixT2T <- as.matrix(T2T_new)
    assign("matrixT2T", matrixT2T)
  }
}

## Save processed matrices ---------------------------------------------------
saveRDS(matrixhg18,   "hg18_for_wilcox.RDS")
saveRDS(matrixhg19,   "hg19_for_wilcox.RDS")
saveRDS(matrixv30,    "v30_for_wilcox.RDS")
saveRDS(matrixv43,    "v43_for_wilcox.RDS")
saveRDS(matrixT2T,    "T2T_for_wilcox.RDS")
