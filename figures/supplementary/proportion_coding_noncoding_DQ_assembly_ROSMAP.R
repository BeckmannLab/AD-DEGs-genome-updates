# ==== Libraries ====
library(ggplot2)

# ==== Load Preprocessed Data ====
all_info_ROSMAP <- readRDS(input_rds)

# ==== Format Variables for Plotting ====
all_info_ROSMAP$for_plot <- factor(
  all_info_ROSMAP$for_plot,
  levels = c("New_TRUE", "New_FALSE", "Old_TRUE", "Old_FALSE", "notSignif_TRUE")
)

cbPalette <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080", "#999999")

# Reformat assembly labels
assembly_map <- c(
  "GRCh37-NCBI36" = "GRCh37/NCBI36",
  "GRCh38.12-NCBI36" = "GRCh38.p12/NCBI36",
  "GRCh38.13-NCBI36" = "GRCh38.p13/NCBI36",
  "CHM13v2.0-NCBI36" = "T2T-CHM13v2.0/NCBI36",
  "GRCh38.12-GRCh37" = "GRCh38.p12/GRCh37",
  "GRCh38.13-GRCh37" = "GRCh38.p13/GRCh37",
  "CHM13v2.0-GRCh37" = "T2T-CHM13v2.0/GRCh37",
  "GRCh38.13-GRCh38.12" = "GRCh38.p13/GRCh38.p12",
  "CHM13v2.0-GRCh38.12" = "T2T-CHM13v2.0/GRCh38.p12",
  "CHM13v2.0-GRCh38.13" = "T2T-CHM13v2.0/GRCh38.p13"
)

all_info_ROSMAP$assembly <- factor(
  recode(all_info_ROSMAP$assembly, !!!assembly_map),
  levels = unname(assembly_map)
)

# ==== Create Plot ====
g <- ggplot(all_info_ROSMAP, aes(x = assembly, fill = for_plot)) +
  geom_bar(stat = "count", position = "dodge") +
  facet_wrap(
    ~status,
    labeller = labeller(status = c(non_coding = "Non-Coding", protein_coding = "Protein-Coding"))
  ) +
  theme_minimal() +
  xlab("Assembly") +
  ylab("Genes Expressed") +
  scale_fill_manual(
    values = cbPalette,
    labels = c(
      "common genes LogFC > 0",
      "unique genes LogFC > 0",
      "common genes LogFC ≤ 0",
      "unique genes LogFC ≤ 0",
      "common genes FDR > 0.05"
    )
  ) +
  labs(fill = "ROSMAP Differential Quantification") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x.bottom = element_text(size = 15),
    axis.title.y.left = element_text(size = 15),
    legend.position = c(0.18, 0.8),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    strip.text = element_text(size = 12)
  )

# ==== Save Plot ====
ggsave(
  filename = output_pdf,
  plot = g,
  device = cairo_pdf,
  width = 12,
  height = 6
)
