# Load libraries
library(corrplot)
library(RColorBrewer)

# Read in data
rosmap <- readRDS(rosmap_path)
msbb   <- readRDS(msbb_path)

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

# Descriptive assembly names
rename_labels <- function(labs) {
  map <- c(
    "hg19_hg18" = "GRCh37/NCBI36",
    "v30_hg18"  = "GRCh38.p12/NCBI36",
    "v43_hg18"  = "GRCh38.p13/NCBI36",
    "v30_hg19"  = "GRCh38.p12/GRCh37",
    "v43_hg19"  = "GRCh38.p13/GRCh37",
    "v43_v30"   = "GRCh38.p13/GRCh38.p12",
    "T2T_hg18"  = "T2T-CHM13v2.0/NCBI36",
    "T2T_hg19"  = "T2T-CHM13v2.0/GRCh37",
    "T2T_v30"   = "T2T-CHM13v2.0/GRCh38.p12",
    "T2T_v43"   = "T2T-CHM13v2.0/GRCh38.p13"
  )
  for (pat in names(map)) {
    labs <- gsub(pat, map[pat], labs)
  }
  labs
}
rownames(res) <- sub("_rosmap$", "", sub("^logFC_", "", rename_labels(rownames(res))))
colnames(res) <- sub("_msbb$",   "", sub("^logFC_", "", rename_labels(colnames(res))))

# Diverging palette
div_col    <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(200)
div_breaks <- seq(-1, 1, length.out = length(div_col) + 1)

# Plot
pdf(output_pdf, width = 7, height = 7)

# **increase right margin** so we have room to shift legend‐labels right
par(mar = c(5, 5, 2, 6))

corrplot(
  res,
  type         = "upper",
  method       = "color",
  col          = div_col,
  breaks       = div_breaks,
  addgrid.col  = "white",
  diag         = TRUE,
  is.corr      = FALSE,
  tl.col       = "black",
  tl.srt       = 45,
  tl.cex       = 0.8,
  number.digits= 2,
  cl.cex       = 0.7,
  cl.length    = 11,    # ticks every 0.2
  cl.offset    = 3      # shift legend (bar + labels) further right
)

# Axis labels
mtext("MSBB logFC",   side = 1, line = 3, cex = 1)
mtext("ROSMAP logFC", side = 2, line = 3, cex = 1)

dev.off()
