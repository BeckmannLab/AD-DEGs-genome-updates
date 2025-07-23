# Load libraries
library(corrplot)
library(RColorBrewer)

# Read in data
rosmap <- readRDS("with_signif_categories_results_assembly_ROSMAP_1.30.24.RDS") #see AD-DEGs-genome-updates/misc/with_signif_categories_results_assembly_ROSMAP_1.30.24.R
msbb   <- readRDS("with_signif_categories_results_assembly_MSBB_3.30.24.RDS") #see AD-DEGs-genome-updates/misc/with_signif_categories_results_assembly_MSBB_3.30.24.R

# Helper to pull out logFC & p-values
process_dataset <- function(df, prefix, select_cols = 1:61) {
  df_sub <- df[, select_cols]
  rownames(df_sub) <- df_sub$gene_id

  logfc_cols <- grep("logFC", colnames(df_sub), value = TRUE)
  logfc      <- na.omit(df_sub[, logfc_cols])
  colnames(logfc) <- paste0(logfc_cols, "_", prefix)

  pval_cols <- grep("adj.P.Val", colnames(df_sub), value = TRUE)
  adjp      <- na.omit(df_sub[, pval_cols])
  colnames(adjp) <- paste0(pval_cols, "_", prefix)

  list(logfc = logfc, adjp = adjp)
}

# Process both datasets
rosmap_res <- process_dataset(rosmap, "rosmap")
msbb_res   <- process_dataset(msbb,   "msbb")

# Merge logFC by gene_id
logfc_merged <- merge(
  msbb_res$logfc,
  rosmap_res$logfc,
  by = "row.names", all = TRUE
)
rownames(logfc_merged) <- logfc_merged$Row.names
logfc_merged$Row.names <- NULL

# Build all pair‐wise combos and compute Spearman ρ
combos <- expand.grid(
  row = colnames(rosmap_res$logfc),
  col = colnames(msbb_res$logfc),
  stringsAsFactors = FALSE
)
res <- matrix(
  NA,
  nrow = ncol(rosmap_res$logfc),
  ncol = ncol(msbb_res$logfc),
  dimnames = list(
    colnames(rosmap_res$logfc),
    colnames(msbb_res$logfc)
  )
)
for (i in seq_len(nrow(combos))) {
  rnm <- combos$row[i]
  cnm <- combos$col[i]
  tmp <- na.omit(logfc_merged[, c(rnm, cnm)])
  res[rnm, cnm] <- cor(tmp[,1], tmp[,2], method = "spearman")
}

saveRDS(res, "cor_msbb_rosmap.RDS")
