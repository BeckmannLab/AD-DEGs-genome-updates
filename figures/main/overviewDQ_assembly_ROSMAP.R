# ==== Set Working Directories (Edit or Remove Before Publishing) ====
data_dir <- "/sc/arion/projects/mscic1/results/anina/fun_project_4.23/analysis/"
output_dir <- "/hpc/users/lunda02/www/plots/assembly_paper/assembly/"

input_rds_all <- file.path(data_dir, "with_signif_categories_results_assembly_ROSMAP_1.30.24.RDS")
input_rds_msbb <- file.path(data_dir, "all_info_MSBB.RDS")
input_rds_rosmap <- file.path(data_dir, "all_info_ROSMAP.RDS")
output_pdf <- file.path(output_dir, "overviewDQ_assembly_ROSMAP.pdf")

# ==== Libraries ====
library(ggplot2)

# ==== Load Data ====
all <- readRDS(input_rds_all)
all_info_MSBB <- readRDS(input_rds_msbb)
info_all_ROSMAP <- readRDS(input_rds_rosmap)

# ==== Set Factor Levels and Colors ====
info_all_ROSMAP$for_plot <- factor(
  info_all_ROSMAP$for_plot,
  levels = c("New_TRUE", "New_FALSE", "Old_TRUE", "Old_FALSE", "notSignif_TRUE")
)

info_all_ROSMAP$assembly <- factor(
  info_all_ROSMAP$assembly,
  levels = c(
    "GRCh37-NCBI36", "GRCh38.12-NCBI36", "GRCh38.13-NCBI36", "CHM13v2.0-NCBI36",
    "GRCh38.12-GRCh37", "GRCh38.13-GRCh37", "CHM13v2.0-GRCh37",
    "GRCh38.13-GRCh38.12", "CHM13v2.0-GRCh38.12", "CHM13v2.0-GRCh38.13"
  )
)

cbPalette <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080", "#999999")

# ==== Create Plot ====
g <- ggplot(info_all_ROSMAP, aes(x = assembly, fill = for_plot)) +
  geom_bar(stat = "count", position = "dodge") +
  theme_minimal() +
  xlab("Assembly") +
  ylab("Genes Expressed") +
  scale_fill_manual(
    values = cbPalette,
    labels = c(
      "common genes LogFC > 0",
      "unique genes LogFC > 0",
      "common genes LogFC \u2264 0",
      "unique genes LogFC \u2264 0",
      "common genes FDR > 0.05"
    )
  ) +
  labs(fill = "ROSMAP Differential Quantification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x.bottom = element_text(size = 22),
    axis.title.y.left = element_text(size = 22),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17)
  )

# ==== Save Plot ====
ggsave(
  filename = output_pdf,
  plot = g,
  device = cairo_pdf,
  width = 10,
  height = 8
)
