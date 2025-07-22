# -----------------------------
# Clear workspace
# -----------------------------
rm(list=ls())  # Remove all existing objects for a clean session

# -----------------------------
# Load libraries
# -----------------------------

library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)

# -----------------------------
# Custom functions
# -----------------------------

# add_decimal_numbers: append .1, .2, ... to duplicate gene IDs
# - target_values: vector of unique IDs to annotate
# - value_list: full list of IDs (including duplicates)
add_decimal_numbers <- function(target_values, value_list) {
    target_indices <- match(value_list, target_values)
    modified_list <- ifelse(
        !is.na(target_indices),
        paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
        value_list
    )
    return(modified_list)
}


# duplicate_rows_based_on_count: replicate rows in dataset based on count_dataset frequencies
# - dataset: data.frame to duplicate rows in
# - column_to_duplicate: column name indicating ID to replicate
# - count_dataset: data.frame with counts of IDs
duplicate_rows_based_on_count <- function(dataset, column_to_duplicate, count_dataset) {
    # Determine overlapping IDs
    overlap <- intersect(dataset[, column_to_duplicate], count_dataset[, column_to_duplicate])
    # Table of counts for overlapping IDs
    count_values <- table(count_dataset$Geneid_noVersion[count_dataset$Geneid_noVersion %in% overlap])
    # Identify IDs not duplicated
    values_not_in_count <- setdiff(dataset$Geneid_noVersion, names(count_values))

    # Subset of genes with single occurrence
    dataset_1gene <- dataset[match(names(which(count_values == 1)), dataset$Geneid_noVersion), ]
    # Subset of genes not in count dataset
    dataset_non_count <- dataset[match(values_not_in_count, dataset$Geneid_noVersion), ]
    # Subset of genes with multiple occurrences
    dataset_more <- dataset[match(names(which(count_values > 1)), dataset$Geneid_noVersion), ]
    assert_that(sum(is.na(dataset_more)) == 0)

    # Duplicate rows per occurrence count
    for (i in seq_len(nrow(dataset_more))) {
        value <- dataset_more[i, column_to_duplicate]
        occurrence_count <- count_values[value]
        dup <- dataset_more[rep(i, occurrence_count), ]
        # Append .1, .2, etc to rownames and ID
        rownames(dup) <- dup$Geneid_noVersion <- paste0(value, ".", seq_len(occurrence_count))
        if (i == 1) {
            dataset_dup <- dup
        } else {
            dataset_dup <- rbind(dataset_dup, dup)
        }
    }
    # Combine all subsets back together
    rbind(dataset_1gene, dataset_dup, dataset_non_count)
}

# -----------------------------
# Load count matrices for each assembly
# -----------------------------
# hg18:
genedata_org_hg18 <- fread("allcount_matrix_2023-04-11.txt", data.table = FALSE)
rownames(genedata_org_hg18) <- genedata_org_hg18$Geneid
genedata_org_hg18$Geneid <- NULL

# hg19:
genedata_org_hg19 <- fread("allcount_matrix_2023-04-12.txt", data.table = FALSE)
rownames(genedata_org_hg19) <- genedata_org_hg19$Geneid
genedata_org_hg19$Geneid <- NULL

# v30:
genedata_org_v30 <- fread("allcount_matrix_2023-04-11.txt", data.table = FALSE)
rownames(genedata_org_v30) <- genedata_org_v30$Geneid
genedata_org_v30$Geneid <- NULL

# v43:
genedata_org_v43 <- fread("allcount_matrix_2023-04-11.txt", data.table = FALSE)
rownames(genedata_org_v43) <- genedata_org_v43$Geneid
genedata_org_v43$Geneid <- NULL

# T2T:
genedata_org_T2T <- fread("allcount_matrix_2023-07-03.txt", data.table = FALSE)
rownames(genedata_org_T2T) <- genedata_org_T2T$Geneid
genedata_org_T2T$Geneid <- NULL

# -----------------------------
# Mapping file for T2T assembly
# -----------------------------
map <- fread("map_genename_T2T_ensembl_from_gtf.txt")

# -----------------------------
# Pairwise mapping grid (non-T2T)
# -----------------------------
v1 <- c("genedata_org_hg18","genedata_org_hg18","genedata_org_v30","genedata_org_hg19","genedata_org_hg19","genedata_org_v30")
v2 <- c("hg18_name","hg18_name","v30_name","hg19_name","hg19_name","v30_name")
v3 <- c("hg18","hg18","v30","hg19","hg19","v30")
v4 <- c("genedata_org_hg19","genedata_org_v30","genedata_org_v43","genedata_org_hg19","genedata_org_v30","genedata_org_v43")
v5 <- c("hg19_name","v30_name","v43_name","hg19_name","v30_name","v43_name")
v6 <- c("hg19","v30","v43","hg19","v30","v43")

grid <- as.data.frame(rbind(v1, v2, v3, v4, v5, v6), stringsAsFactors = FALSE)

# Pre-allocate list for results
results_list <- vector("list", ncol(grid))

# Loop over non-T2T pairs
for (i in seq_len(ncol(grid))) {
    # --- Assembly 1 reference
    ds1 <- get(grid[1, i])  # Load dataset by name
    gd1 <- ds1[, 6:ncol(ds1)]  # Keep only count columns
    # Strip version suffix from rownames
    gd1$Geneid_noVersion <- sapply(strsplit(rownames(gd1), "\\.", fixed = FALSE), function(x) {
        if (!grepl("_", tail(x, 1))) {
            x[1]
        } else {
            parts <- unlist(strsplit(x[2], "_"))
            paste(x[1], paste(parts[2:3], collapse = "_"), sep = "_")
        }
    })
    var1 <- grid[2, i]  # Column name for original IDs
    gd1[[var1]] <- rownames(gd1)
    ref1 <- gd1[, c(var1, "Geneid_noVersion")]
    colnames(ref1) <- c("assembly_1_original_gene_name", "common_name")
    ref1$assembly_1_name <- grid[3, i]

    # --- Assembly 2 reference
    ds2 <- get(grid[4, i])
    gd2 <- ds2[, 6:ncol(ds2)]
    # Strip version suffix
    gd2$Geneid_noVersion <- sapply(strsplit(rownames(gd2), "\\.", fixed = FALSE), function(x) {
        if (!grepl("_", tail(x, 1))) {
            x[1]
        } else {
            parts <- unlist(strsplit(x[2], "_"))
            paste(x[1], paste(parts[2:3], collapse = "_"), sep = "_")
        }
    })
    var2 <- grid[5, i]
    gd2[[var2]] <- rownames(gd2)
    ref2 <- gd2[, c(var2, "Geneid_noVersion")]
    colnames(ref2) <- c("assembly_2_original_gene_name", "common_name")
    ref2$assembly_2_name <- grid[6, i]

    # --- Merge references by common gene name
    both <- merge(ref1, ref2, by = "common_name", all = TRUE)
    both$assembly_1_name <- grid[3, i]
    both$assembly_2_name <- grid[6, i]
    results_list[[i]] <- both
}

# Combine all non-T2T results and add combo column
final_results <- do.call(rbind, results_list)
final_results$combo <- paste(final_results$assembly_1_name, final_results$assembly_2_name, sep = "_")

# -----------------------------
# Pairwise mapping including T2T
# -----------------------------
# Define grid for T2T comparisons
v1_t2t <- c("genedata_org_hg18", "genedata_org_hg19", "genedata_org_v30", "genedata_org_v43")
v2_t2t <- c("hg18_name", "hg19_name", "v30_name", "v43_name")
v3_t2t <- c("hg18", "hg19", "v30", "v43")

grid2 <- as.data.frame(rbind(v1_t2t, v2_t2t, v3_t2t), stringsAsFactors = FALSE)
# Pre-allocate list
results_list2 <- vector("list", ncol(grid2))

for (i in seq_len(ncol(grid2))) {
    print(i)
    # Load T2T data and attach gene IDs
    genedata <- genedata_org_T2T[, 6:ncol(genedata_org_T2T)]
    genedata$gene_id <- rownames(genedata)
    assert_that(nrow(genedata) == nrow(map))  # Ensure same number of rows
    genedata <- merge(genedata, map, by.x = "gene_id", by.y = "gene_id")

    # Create common_name by stripping version and merging parent ID
    genedata$Geneid_noVersion <- unlist(lapply(strsplit(as.character(genedata$projection_parent_gene), ".", fixed = TRUE), function(x) {
        if (sum(grepl("_", x)) == FALSE) {
            x[1]
        } else {
            paste(x[1], paste(unlist(strsplit(as.character(x[2]), "_"))[2:3], collapse = "_"), sep = "_")
        }
    }))

    # Identify duplicated parent genes
    dupl <- genedata[duplicated(genedata$projection_parent_gene), ]
    uni_ids <- unique(dupl$projection_parent_gene)
    uni_ids <- data.frame(ids = uni_ids, stringsAsFactors = FALSE)
    uni_ids$Geneid_noVersion <- unlist(lapply(strsplit(as.character(uni_ids$ids), ".", fixed = TRUE), function(x) x[1]))

    # Append decimals to duplicates
    values <- genedata$Geneid_noVersion
    targets <- uni_ids$Geneid_noVersion
    genedata$modified_values <- add_decimal_numbers(targets, values)

    # Set rownames to modified values
    rownames(genedata) <- genedata$modified_values
    genedata_formod <- genedata

    # Prepare T2T reference
    T2T_ref <- genedata[, c("gene_id", "modified_values")]
    colnames(T2T_ref) <- c("assembly_2_original_gene_name", "common_name")
    T2T_ref$assembly_2_name <- "T2T"

    # --- Compare with other assemblies
    ds1 <- get(grid2[1, i])  # Assembly dataset
    gd1 <- ds1[, 6:ncol(ds1)]
    # Strip version
    gd1$Geneid_noVersion <- sapply(strsplit(rownames(gd1), "\\.", fixed = FALSE), function(x) {
        if (!grepl("_", tail(x, 1))) {
            x[1]
        } else {
            parts <- unlist(strsplit(x[2], "_"))
            paste(x[1], paste(parts[2:3], collapse = "_"), sep = "_")
        }
    })
    var1 <- grid2[2, i]
    gd1[[var1]] <- rownames(gd1)

    # Duplicate rows to match T2T counts
    duplicated_df <- duplicate_rows_based_on_count(gd1, "Geneid_noVersion", genedata_formod)

    # Build reference for assembly 1
    reference_all <- data.frame(
        assembly_1_original_gene_name = duplicated_df[[var1]],
        common_name = duplicated_df$Geneid_noVersion,
        stringsAsFactors = FALSE
    )
    reference_all$assembly_1_name <- grid2[3, i]

    # Merge with T2T reference
    reference2_T2T <- merge(reference_all, T2T_ref, by = "common_name", all = TRUE)
    reference2_T2T$assembly_2_name <- "T2T"
    results_list2[[i]] <- reference2_T2T
}

# Combine all T2T results and add combo column
final_results2 <- do.call(rbind, results_list2)
final_results2$combo <- paste(final_results2$assembly_1_name, final_results2$assembly_2_name, sep = "_")

# -----------------------------
# Combine all mapping results
# -----------------------------
all <- rbind(final_results, final_results2)

saveRDS(all,"rbind_map_between_assemblies_10.04.23.RDS")

