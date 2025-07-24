
rm(list=ls()) 
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)

##function to add decimal 
add_decimal_numbers <- function(target_values, value_list) {
  # Find indices of target values in the value list
  target_indices <- match(value_list, target_values)
  
  # Generate modified values for target values
  modified_list <- ifelse(
    !is.na(target_indices),
    paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
    value_list
  )
  
  return(modified_list)
}


#########################
#wilcox plots
#########################

genedata_org_hg18=fread("allcount_matrix_2023-04-11.txt", data.table=FALSE) #counts
rownames(genedata_org_hg18)=genedata_org_hg18$Geneid
genedata_org_hg18$Geneid=NULL
dim(genedata_org_hg18)

genedata_org_hg19=fread("allcount_matrix_2023-04-12.txt",data.table=FALSE) #counts
rownames(genedata_org_hg19)=genedata_org_hg19$Geneid
genedata_org_hg19$Geneid=NULL
dim(genedata_org_hg19)

genedata_org_v30=fread("allcount_matrix_2023-04-11.txt", data.table=FALSE) #counts
rownames(genedata_org_v30)=genedata_org_v30$Geneid
genedata_org_v30$Geneid=NULL
dim(genedata_org_v30)

genedata_org_v43=fread("allcount_matrix_2023-04-11.txt", data.table= FALSE) #counts
rownames(genedata_org_v43)=genedata_org_v43$Geneid
genedata_org_v43$Geneid=NULL
dim(genedata_org_v43)

genedata_org_T2T=fread("allcount_matrix_2023-07-03.txt", data.table = FALSE) #counts

map <- fread("map_genename_T2T_ensembl_from_gtf.txt") #see AD-DEGs-updates/files/map_genename_T2T_ensembl_from_gtf.txt

assembly_map <- readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-updates/misc/rbind_map_between_assemblies_10.04.23.R


v1 = c("genedata_org_hg18","genedata_org_hg19","genedata_org_v30","genedata_org_v43", "genedata_org_T2T")
v2 = c("hg18","hg19", "v30","v43", "T2T")
grid = as.data.frame(cbind(v1, v2))


for (x in seq_len(nrow(grid))){
  print(grid[x,2])
  if (grid[x,2] != "T2T"){
    print(grid[x,2])
    og_ref = get(grid[x,1])
    ref=og_ref[,6:644]
    ref_isexpr = rowSums(cpm(ref)>=1) >= 0.1*ncol(ref) 
    ref=ref[ref_isexpr,]
    ref$Geneid_noVersion=unlist(lapply(strsplit(as.character(rownames(ref)),".",fixed=TRUE),function(x){
          if(sum(grepl("_", x))==FALSE){
            x[1]
          }else{
            paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
          }
          }))
    rownames(ref)=ref$Geneid_noVersion
    ref$Geneid_noVersion=NULL
    ref <- as.matrix(ref)
    assign(paste0("matrix",grid[x,2]),ref)
}else{
  og_ref = get(grid[x,1])
  T2T <- merge(og_ref, map, by.x = "Geneid" , by.y = "gene_id")
  T2T$Geneid_noVersion=unlist(lapply(strsplit(as.character(T2T$projection_parent_gene),".",fixed=TRUE),function(x){
        if(sum(grepl("_", x))==FALSE){
          x[1]
        }else{
          paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
        }
        }))
  ##see which genes are duplicated 
    dupl <- T2T[duplicated(T2T$projection_parent_gene),]

    #unique ids
    uni_ids <- unique(dupl$projection_parent_gene)
    uni_ids <- as.data.frame(uni_ids)
    colnames(uni_ids) <- "ids"
    uni_ids$Geneid_noVersion=unlist(lapply(strsplit(as.character(uni_ids$ids),".",fixed=TRUE),function(x){x[1]}))

    #add decimal 
    values <- T2T$Geneid_noVersion  # Vector of values to modify
    targets <- uni_ids$Geneid_noVersion  # Vector of target values

    modified_values <- add_decimal_numbers(targets, values)
    T2T$modified_values <- unlist(modified_values)

    #make rownames 
    rownames(T2T) <- T2T$modified_values
    T2T$projection_parent_gene=NULL
    T2T$Geneid_noVersion=NULL
    T2T$modified_values=NULL

    T2T_new=T2T[,7:645]

    T2T_new_isexpr = rowSums(cpm(T2T_new)>=1) >= 0.1*ncol(T2T_new) 
    T2T_new=T2T_new[T2T_new_isexpr,]

    matrixT2T <- as.matrix(T2T_new)
    assign("matrixT2T",matrixT2T)
}}

saveRDS(matrixhg18, "hg18_for_wilcox.RDS")
saveRDS(matrixhg19, "hg19_for_wilcox.RDS")
saveRDS(matrixv30, "v30_for_wilcox.RDS")
saveRDS(matrixv43, "v43_for_wilcox.RDS")
saveRDS(matrixT2T, "T2T_for_wilcox.RDS")
