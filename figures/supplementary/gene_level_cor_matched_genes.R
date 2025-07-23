##LIBRARIES
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(reshape2)
  library(variancePartition)
  library(assertthat)
  library(gridExtra)
  library(tidyverse)
})

##LOAD
path_hg18 <- "allcount_matrix_2023-04-11.txt"
path_hg19 <- "allcount_matrix_2023-04-12.txt"
path_v30 <- "allcount_matrix_2023-04-11.txt"
path_v43 <- "allcount_matrix_2023-04-11.txt"
path_T2T <- "allcount_matrix_2023-07-03.txt"
path_map <- "map_genename_T2T_ensembl_from_gtf.txt" # see AD_DEGs-genome-updates/files/map_genename_T2T_ensembl_from_gtf.txt
path_assembly_map <- "rbind_map_between_assemblies_10.04.23.RDS" #see AD_DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.RDS


##FUNCTIONS
load_expression <- function(path, gene_col = "Geneid", data_cols = 6:644) {
  df <- fread(path, data.table = FALSE)
  rownames(df) <- df[[gene_col]]
  df[[gene_col]] <- NULL
  df[, data_cols]
}

filter_expressed_genes <- function(expr_matrix, min_cpm = 1, min_prop = 0.1) {
  keep <- rowSums(cpm(expr_matrix) >= min_cpm) >= (min_prop * ncol(expr_matrix))
  expr_matrix[keep, ]
}

strip_gene_version <- function(genes) {
  sapply(strsplit(genes, ".", fixed = TRUE), function(x) {
    if (!any(grepl("_", x))) x[1] else paste(x[1], paste(unlist(strsplit(x[2], "_"))[2:3], collapse = "_"), sep = "_")
  })
}

map_T2T_to_assembly <- function(genedata_T2T, assembly_map, combo_label) {
  mapping <- assembly_map[assembly_map$combo == combo_label, ]
  merged <- merge(genedata_T2T, mapping, by.x = "Geneid", by.y = "assembly_2_orginal_gene_name")
  rownames(merged) <- merged$common_name
  keep_cols <- setdiff(colnames(merged), c("common_name", "assembly_1_orginal_gene_name", "assembly_1_name", "combo", "assembly_2_name"))
  expressed <- merged[, keep_cols]
  filter_expressed_genes(expressed[, 7:645])
}

compute_correlation_plot <- function(df1, df2, label, title) {
  common_genes <- intersect(rownames(df1), rownames(df2))
  corrs <- mapply(cor, as.data.frame(t(df1[common_genes, ])), as.data.frame(t(df2[common_genes, ])))
  corr_df <- melt(as.data.frame(corrs))
  corr_df$assembly <- label
  gg <- ggplot(corr_df, aes(x = value)) +
    geom_density(fill = "skyblue", color = "darkblue", alpha = 0.7) +
    labs(title = title, x = "Correlation Values", y = "Density") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 10))
  list(data = corr_df, plot = gg)
}

##LOAD DATA
hg18 <- load_expression(path_hg18)
hg19 <- load_expression(path_hg19)
v30 <- load_expression(path_v30)
v43 <- load_expression(path_v43)
genedata_T2T <- fread(path_T2T, data.table = FALSE)
rownames(genedata_T2T) <- genedata_T2T$Geneid

assembly_map <- readRDS(path_assembly_map)

##FILTER FOR EXPRESSED GENES
hg18 <- filter_expressed_genes(hg18)
hg19 <- filter_expressed_genes(hg19)
v30 <- filter_expressed_genes(v30)
v43 <- filter_expressed_genes(v43)

##STRIP VERSIONS FROM IDs
rownames(hg18) <- strip_gene_version(rownames(hg18))
rownames(hg19) <- strip_gene_version(rownames(hg19))
rownames(v30) <- strip_gene_version(rownames(v30))
rownames(v43) <- strip_gene_version(rownames(v43))

##MAP T2T TO OTHER ASSEMBLIES
T2T_hg18 <- map_T2T_to_assembly(genedata_T2T, assembly_map, "hg18_T2T")
T2T_hg19 <- map_T2T_to_assembly(genedata_T2T, assembly_map, "hg19_T2T")
T2T_v30 <- map_T2T_to_assembly(genedata_T2T, assembly_map, "v30_T2T")
T2T_v43 <- map_T2T_to_assembly(genedata_T2T, assembly_map, "v43_T2T")

##COMPUTE CORRELATIONS
pairs <- list(
  list(hg18, hg19, "hg19_hg18", "Correlation: GRCh37 and NCBI36"),
  list(hg18, v30, "v30_hg18", "Correlation: GRCh38.p12 and NCBI36"),
  list(hg18, v43, "v43_hg18", "Correlation: GRCh38.p13 and NCBI36"),
  list(hg18, T2T_hg18, "T2T_hg18", "Correlation: CHM13v2.0 and NCBI36"),
  list(hg19, v30, "v30_hg19", "Correlation: GRCh38.p12 and GRCh37"),
  list(hg19, v43, "v43_hg19", "Correlation: GRCh38.p13 and GRCh37"),
  list(hg19, T2T_hg19, "T2T_hg19", "Correlation: CHM13v2.0 and GRCh37"),
  list(v30, v43, "v43_v30", "Correlation: GRCh38.p13 and GRCh38.p12"),
  list(v30, T2T_v30, "T2T_v30", "Correlation: CHM13v2.0 and GRCh38.p12"),
  list(v43, T2T_v43, "T2T_v43", "Correlation: CHM13v2.0 and GRCh38.p13")
)

results <- lapply(pairs, function(p) compute_correlation_plot(p[[1]], p[[2]], p[[3]], p[[4]]))

##COMBINE PLOT DATA
all_corrs <- do.call(rbind, lapply(results, `[[`, "data"))

##RENAME FOR PLOT LABELS
assembly_labels <- c(
  "hg19_hg18" = "GRCh37/NCBI36",
  "v30_hg18" = "GRCh38.p12/NCBI36",
  "v43_hg18" = "GRCh38.p13/NCBI36",
  "T2T_hg18" = "T2T-CHM13v2.0/NCBI36",
  "v30_hg19" = "GRCh38.p12/GRCh37",
  "v43_hg19" = "GRCh38.p13/GRCh37",
  "T2T_hg19" = "T2T-CHM13v2.0/GRCh37",
  "v43_v30" = "GRCh38.p13/GRCh38.p12",
  "T2T_v30" = "T2T-CHM13v2.0/GRCh38.p12",
  "T2T_v43" = "T2T-CHM13v2.0/GRCh38.p13"
)

all_corrs$assembly <- factor(assembly_labels[all_corrs$assembly],
  levels = assembly_labels)

##PLOT ALL
g <- ggplot(all_corrs, aes(x = value, color = assembly)) +
  geom_density(aes(y = ..scaled..), alpha = 0.3) +
  labs(
       x = "Correlation Values", y = "Density", color = "Reference") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 15)) +
  xlim(0.999, 1)

ggsave(g, filename = plot_outfile, useDingbats = FALSE)

