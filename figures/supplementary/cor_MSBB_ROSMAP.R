#Supplementary figure 5B
################################################################################
# LOAD LIBRARIES
################################################################################
library(ggplot2)

################################################################################
# DEFINE FUNCTIONS
################################################################################
# Convert assembly comparison codes into human-readable labels
rename_labels <- function(labels) {
  mapping <- c(
    hg19_hg18 = "GRCh37-NCBI36",
    v30_hg18  = "GRCh38.p12-NCBI36",
    v43_hg18  = "GRCh38.p13-NCBI36",
    v30_hg19  = "GRCh38.p12-GRCh37",
    v43_hg19  = "GRCh38.p13-GRCh37",
    v43_v30   = "GRCh38.p13-GRCh38.p12",
    T2T_hg18  = "T2T-CHM13v2.0-NCBI36",
    T2T_hg19  = "T2T-CHM13v2.0-GRCh37",
    T2T_v30   = "T2T-CHM13v2.0-GRCh38.p12",
    T2T_v43   = "T2T-CHM13v2.0-GRCh38.p13"
  )
  sapply(labels, function(x) {
    if (x %in% names(mapping)) mapping[[x]] else x
  }, USE.NAMES = FALSE)
}

################################################################################
# LOAD DATA
################################################################################
res <- readRDS("cor_msbb_rosmap.RDS") #see AD-DEGs-genome-updates/misc/cor_msbb_rosmap.R

################################################################################
# PREPARE DIAGONAL DATA FRAME
################################################################################
# Extract diagonal entries (matching indices) for ROSMAP vs MSBB correlations
i   <- seq_len(min(nrow(res), ncol(res)))
diag_df <- data.frame(
  rosmaps = rownames(res)[i],
  msbb    = colnames(res)[i],
  corr    = res[cbind(i, i)],
  stringsAsFactors = FALSE
)

# Strip metric suffix and 'logFC_' prefix, then map to readable labels
diag_df$reference <- sub("^logFC_", "", sub("_[^_]+$", "", diag_df$rosmaps))
diag_df$reference <- rename_labels(diag_df$reference)
diag_df$metric    <- "LogFC"

# Set factor levels for consistent ordering in plots
level_order <- c(
  "GRCh37-NCBI36",
  "GRCh38.p12-NCBI36",
  "GRCh38.p13-NCBI36",
  "T2T-CHM13v2.0-NCBI36",
  "GRCh38.p12-GRCh37",
  "GRCh38.p13-GRCh37",
  "T2T-CHM13v2.0-GRCh37",
  "GRCh38.p13-GRCh38.p12",
  "T2T-CHM13v2.0-GRCh38.p12",
  "T2T-CHM13v2.0-GRCh38.p13"
)
diag_df$reference <- factor(diag_df$reference, levels = level_order)

################################################################################
# CREATE SCATTER PLOT
################################################################################
g <- ggplot(diag_df, aes(x = reference, y = corr)) +
  geom_point(size = 6) +
  scale_y_continuous(
    limits = c(min(diag_df$corr) - 0.05, max(diag_df$corr) + 0.05)
  ) +
  labs(
    x = "Reference comparison",
    y = "LogFC correlation\n(ROSMAP vs MSBB)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title  = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18)
  )

################################################################################
# SAVE PLOT
################################################################################
ggsave(
  filename    = output_plot,
  plot        = g,
  useDingbats = FALSE,
  width       = 9,
  height      = 7
)

