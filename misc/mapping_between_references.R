## Load Libraries
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(assertthat)
  library(tidyverse); library(BiocParallel)
})

##-----------------------------------------------
## Utility Functions
##-----------------------------------------------

strip_gene_version <- function(ids) {
  sapply(strsplit(as.character(ids), ".", fixed = TRUE), function(x) {
    if (!any(grepl("_", x))) x[1]
    else paste(x[1], paste(unlist(strsplit(x[2], "_"))[2:3], collapse = "_"), sep = "_")
  })
}

add_decimal_numbers <- function(target_values, value_list) {
  target_indices <- match(value_list, target_values)
  ifelse(!is.na(target_indices),
         paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
         value_list)
}

extract_gene_metadata <- function(mat, prefix) {
  gene_data <- mat[, 6:ncol(mat)]
  gene_data <- gene_data %>%
    mutate(Geneid_noVersion = strip_gene_version(rownames(gene_data))) %>%
    mutate(gene_name = rownames(gene_data)) %>%
    as.data.frame()
  setNames(gene_data[, c("gene_name", "Geneid_noVersion")],
           c(paste0(prefix, "_orginal_gene_name"), "common_name"))
}

duplicate_rows_based_on_count <- function(dataset, column_to_duplicate, count_dataset) {
  overlap <- intersect(dataset[[column_to_duplicate]], count_dataset[[column_to_duplicate]])
  count_values <- table(count_dataset[[column_to_duplicate]][count_dataset[[column_to_duplicate]] %in% overlap])
  values_not_in_count <- setdiff(dataset[[column_to_duplicate]], names(count_values))
  
  dataset_1gene <- dataset[match(names(which(count_values == 1)), dataset[[column_to_duplicate]]), ]
  dataset_non_count <- dataset[match(values_not_in_count, dataset[[column_to_duplicate]]), ]
  
  dataset_more <- dataset[match(names(which(count_values > 1)), dataset[[column_to_duplicate]]), ]
  assert_that(sum(is.na(dataset_more)) == 0)
  
  dataset_dup <- do.call(rbind, lapply(seq_len(nrow(dataset_more)), function(i) {
    value <- dataset_more[i, column_to_duplicate]
    occurrence_count <- count_values[value]
    dup <- dataset_more[rep(value, occurrence_count), ]
    rownames(dup) <- dup[[column_to_duplicate]] <- paste0(value, ".", seq_len(occurrence_count))
    dup
  }))
  
  rbind(dataset_1gene, dataset_dup, dataset_non_count)
}

merge_pair <- function(df1, df2, name1, name2) {
  df_merged <- merge(df1, df2, by = "common_name", all = TRUE)
  df_merged$assembly_1_name <- name1
  df_merged$assembly_2_name <- name2
  df_merged
}

##-----------------------------------------------
## Load Data
##-----------------------------------------------

genedata_list <- lapply(expr_files, function(path) {
  dat <- fread(path, data.table = FALSE)
  rownames(dat) <- dat$Geneid
  dat$Geneid <- NULL
  dat
})
map <- fread(map_file)

##-----------------------------------------------
## Generate Pairwise Mappings
##-----------------------------------------------

meta <- list()
for (name in names(genedata_list)[names(genedata_list) != "T2T"]) {
  meta[[name]] <- extract_gene_metadata(genedata_list[[name]], name)
  meta[[name]]$assembly_name <- name
}

## hg18/hg19/v30/v43 pairwise merges
pairs <- list(
  c("hg18", "hg19"), c("hg18", "v30"), c("hg18", "v43"),
  c("hg19", "v30"), c("hg19", "v43"), c("v30", "v43")
)
pairwise_maps <- lapply(pairs, function(p) merge_pair(meta[[p[1]]], meta[[p[2]]], p[1], p[2]))

##-----------------------------------------------
## Handle T2T Special Case
##-----------------------------------------------

# T2T gene map and version stripping
T2T <- genedata_list[["T2T"]][, 6:ncol(genedata_list[["T2T"]])]
T2T$gene_id <- rownames(T2T)
assert_that(nrow(T2T) == nrow(map))

T2T <- merge(T2T, map, by.x = "gene_id", by.y = "gene_id")
T2T$Geneid_noVersion <- strip_gene_version(T2T$projection_parent_gene)
dups <- T2T[duplicated(T2T$projection_parent_gene), ]
targets <- unique(dups$projection_parent_gene)
target_ids <- strip_gene_version(targets)

T2T$modified_values <- add_decimal_numbers(target_ids, T2T$Geneid_noVersion)
rownames(T2T) <- T2T$modified_values
genedata_formod <- T2T
T2T_map <- T2T[, c("gene_id", "modified_values")]
colnames(T2T_map) <- c("assembly_2_orginal_gene_name", "common_name")
T2T_map$assembly_2_name <- "T2T"

## Duplicated matching for T2T
dup_mappings <- lapply(names(expr_files)[names(expr_files) != "T2T"], function(name) {
  dat <- genedata_list[[name]][, 6:ncol(genedata_list[[name]])]
  dat$Geneid_noVersion <- strip_gene_version(rownames(dat))
  dat[[paste0(name, "_name")]] <- rownames(dat)
  dat <- as.data.frame(dat)
  
  duped <- duplicate_rows_based_on_count(dat, "Geneid_noVersion", genedata_formod)
  orig <- duped[[paste0(name, "_name")]]
  dup <- duped$Geneid_noVersion
  df <- data.frame(assembly_1_orginal_gene_name = orig, common_name = dup)
  df$assembly_1_name <- name
  
  merge(df, T2T_map, by = "common_name", all = TRUE)
})

##-----------------------------------------------
## Final Combine and Save
##-----------------------------------------------

all_mappings <- do.call(rbind, c(pairwise_maps, dup_mappings))
all_mappings$combo <- paste(all_mappings$assembly_1_name, all_mappings$assembly_2_name, sep = "_")

saveRDS(all_mappings, output_file)
