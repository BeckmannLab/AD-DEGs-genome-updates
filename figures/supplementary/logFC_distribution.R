# === Libraries ===
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# === Paths ===
data_path <- "with_signif_categories_results_assembly_MSBB_3.30.24.RDS"
gene_map_path <- "gene_ids_ensembl2symbol_fromHUGO_10JUN2020.csv" #see AD-DEGs-genome-updates/files/gene_ids_ensembl2symbol_fromHUGO_10JUN2020.csv

# === Load data and assign label ===
load_data <- function(subpath, label) {
  df <- read.delim(file.path(base_dir, subpath, "DE", "DE_assembly_corrected_for_id_assembly_test.txt"))
  df$assembly <- label
  return(df)
}

# === Define comparisons and labels ===
comparisons <- list(
  hg19_hg18  = "GRCh37-NCBI36",
  v30_hg18   = "GRCh38.12-NCBI36",
  v43_hg18   = "GRCh38.13-NCBI36",
  T2T_hg18   = "CHM13v2.0-NCBI36",
  v30_hg19   = "GRCh38.12-GRCh37",
  v43_hg19   = "GRCh38.13-GRCh37",
  T2T_hg19   = "CHM13v2.0-GRCh37",
  v30_v43    = "GRCh38.13-GRCh38.12",
  T2T_v30    = "CHM13v2.0-GRCh38.12",
  T2T_v43    = "CHM13v2.0-GRCh38.13"
)

# === Load all datasets into a list ===
data_list <- lapply(names(comparisons), function(key) load_data(key, comparisons[[key]]))
names(data_list) <- names(comparisons)

# === Combine into one dataframe ===
all_data <- do.call(rbind, data_list)
all_data$assembly <- factor(all_data$assembly, levels = comparisons)

# === Global color palette ===
my_colors <- brewer.pal(10, "Paired")
names(my_colors) <- levels(all_data$assembly)

# === Plot theme ===
theme_custom <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

# === Plot: Combined density plot ===
all_plot <- ggplot(all_data, aes(x = logFC, color = assembly)) +
  geom_density(aes(y = ..scaled..)) +
  scale_color_manual(values = my_colors) +
  labs(
    title = "Correlation Values for assembly combos gene counts",
    x = "Correlation Values", y = "Density", fill = "Reference Combinations"
  ) +
  theme_custom +
  theme(
    legend.position = c(0.13, 0.8),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 12),
    axis.title.x.bottom = element_text(size = 15),
    axis.title.y.left = element_text(size = 15)
  )

# === Individual density plot ===
make_density_plot <- function(df, title) {
  ggplot(df, aes(x = logFC)) +
    geom_density(fill = "skyblue", color = "darkblue", alpha = 0.7) +
    labs(title = title, x = "LogFC", y = "Density") +
    theme_custom
}

# === Create all individual plots ===
plot_list <- mapply(
  make_density_plot,
  data_list,
  paste("LogFC", comparisons),
  SIMPLIFY = FALSE
)

# === Save to PDF ===
pdf(output_pdf)
do.call(grid.arrange, c(plot_list, ncol = 2, nrow = 5))
dev.off()
