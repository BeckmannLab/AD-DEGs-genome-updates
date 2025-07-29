#Supplementary Figure 18

rm(list = ls())
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Load data
rosmap <- readRDS("with_signif_categories_results_interaction_ROSMAP_4.10.24.RDS") #see AD-DEGs-genome-updates/misc/with_signif_categories_results_interaction_ROSMAP_4.10.24.R
msbb   <- readRDS("with_signif_categories_results_interaction_MSBB_4.7.24.RDS") #see AD-DEGs-genome-updates/misc/with_signif_categories_results_interaction_MSBB_4.7.24.R

# Remove unwanted columns
exclude_patterns <- list(
  rosmap = c("cogdx_defvsctl", "ceradsc_simplified"),
  msbb   = c("NP.1_test.txt", "CDR_defvsctl", "CERJ_simplified", "PlaqueMean_defvsctl",
             "NP.1_simplified", "PATH.Dx_simplified", "bbscore_simplified", "PATH.Dx")
)

clean_dataset <- function(df, patterns) {
  for (p in patterns) {
    df <- df[, !grepl(p, colnames(df))]
  }
  rownames(df) <- df$gene_id
  return(df)
}

rosmap_all <- clean_dataset(rosmap, exclude_patterns$rosmap)
msbb_all   <- clean_dataset(msbb, exclude_patterns$msbb)

# Standardized column name cleaner
pretty_names <- function(x) {
  x <- gsub("logFC_", "", x)
  x <- gsub("adj.P.Val_", "", x)
  x <- gsub("_[^_]+_[^_]+$", "", x)  # remove _assembly_assembly suffix
  x <- gsub("_", " ", x)
  x <- gsub("simplified", "", x)
  x <- gsub("test.txt", "", x)
  trimws(x)
}

# Clean and extract logFC or adj.P.Val matrix by suffix
clean_names <- function(mat, suffix) {
  mat <- mat[, grepl(suffix, colnames(mat))]
  colnames(mat) <- gsub("\\||-", "_", colnames(mat))
  mat <- na.omit(mat)
  return(mat)
}
# Core correlation analysis function
analyze_correlation_gg <- function(rosmap_df, msbb_df, suffix, label) {
  logfc_ros <- clean_names(rosmap_df[, grep("logFC", colnames(rosmap_df))], suffix)
  logfc_msbb <- clean_names(msbb_df[, grep("logFC", colnames(msbb_df))], suffix)
  
  pval_ros <- clean_names(rosmap_df[, grep("adj.P.Val", colnames(rosmap_df))], suffix)
  pval_msbb <- clean_names(msbb_df[, grep("adj.P.Val", colnames(msbb_df))], suffix)

  both <- merge(logfc_ros, logfc_msbb, by = "row.names", all = TRUE)
  both2 <- merge(pval_ros, pval_msbb, by = "row.names", all = TRUE)
  rownames(both) <- both$Row.names
  rownames(both2) <- both2$Row.names
  both$Row.names <- both2$Row.names <- NULL
  
  comb <- expand.grid(colnames(logfc_ros), colnames(logfc_msbb))
  cor_list <- list()

  for (i in 1:nrow(comb)) {
    e1 <- comb[i, 1]; e2 <- comb[i, 2]
    p1 <- gsub("logFC", "adj.P.Val", e1); p2 <- gsub("logFC", "adj.P.Val", e2)
    
    df <- both[, c(e1, e2)]
    df <- df[complete.cases(df), ]
    df2 <- both2[rownames(df), c(p1, p2)]
    
    df <- df[df2[, 1] <= 1 | df2[, 2] <= 1, ]
    
    if (nrow(df) > 2) {
      cor_val <- cor(df[, 1], df[, 2], method = "spearman")
    } else {
      cor_val <- NA
    }
    
    cor_list[[i]] <- data.frame(
      x = pretty_names(e2),
      y = pretty_names(e1),
      cor = cor_val,
      pair = label
    )
  }
  do.call(rbind, cor_list)
}

# Assembly pairs and labels
assembly_pairs <- c("hg19_hg18", "v30_hg18", "v43_hg18", "T2T_hg18", 
                    "v30_hg19", "v43_hg19", "T2T_hg19", "v43_v30", 
                    "T2T_v30", "T2T_v43")

label_names <- c(
  "hg19_hg18" = "NCBI36/GRCh37",
  "v30_hg18"  = "NCBI36/GRCh38.p12",
  "v43_hg18"  = "NCBI36/GRCh38.p13",
  "T2T_hg18"  = "NCBI36/CHM13v2.0",
  "v30_hg19"  = "GRCh37/GRCh38.p12",
  "v43_hg19"  = "GRCh37/GRCh38.p13",
  "T2T_hg19"  = "GRCh37/CHM13v2.0",
  "v43_v30"   = "GRCh38.p12/GRCh38.p13",
  "T2T_v30"   = "GRCh38.p12/CHM13v2.0",
  "T2T_v43"   = "GRCh38.p13/CHM13v2.0"
)

# Combine results
cor_df_all <- bind_rows(lapply(assembly_pairs, function(pair) {
  analyze_correlation_gg(rosmap_all, msbb_all, pair, label_names[pair])
}))

cor_df_all$y = gsub("braaksc", "Braak",cor_df_all$y )
cor_df_all$y = gsub("ceradsc defvsctl", "CERAD",cor_df_all$y )
cor_df_all$y = gsub("cogdx", "Cognitive",cor_df_all$y )
cor_df_all$y = gsub("cts mmse30", "MMSE30",cor_df_all$y )

cor_df_all$x = gsub("CERJ defvsctl", "CERJ",cor_df_all$x)
cor_df_all$x = gsub("PlaqueMean", "Plaque mean\ndensity",cor_df_all$x)

# Plot: Single figure with facets
plot <- ggplot(cor_df_all, aes(x = x, y = y, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 3, color = "black", na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, na.value = "grey90", limits = c(-1, 1)) +
  labs(fill = expression(atop("Spearman's", rho))) + 
  facet_wrap(~pair, scales = "free", ncol = 3) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(size = 12),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.title = element_blank()
  )

# Save PDF
ggsave(
  filename = "corr_interaction_msbb_rosmap_all_pairs.pdf",
  plot = plot,
  width = 14, height = 10
)
