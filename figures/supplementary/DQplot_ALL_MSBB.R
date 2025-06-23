# Load libraries necessary for data manipulation and plotting
dependency_packages <- c("ggplot2", "data.table", "tidyr", "GGally")
invisible(lapply(dependency_packages, library, character.only = TRUE))

# Read the main results and annotation data into R objects
all           <- readRDS(rds_results)
all_info_MSBB <- readRDS(rds_info)

# Convert the 'for_plot' column into a factor to control the display order in the plot
all_info_MSBB$for_plot <- factor(
  all_info_MSBB$for_plot,
  levels = c("New_TRUE", "New_FALSE", "Old_TRUE", "Old_FALSE", "notSignif_TRUE")
)

# Custom color palette for each differential category passed to scale_fill_manual
cbPalette <- c(
  "#F5793A",  # common genes with positive LogFC
  "#A95AA1",  # unique genes with positive LogFC
  "#85C0F9",  # common genes with non-positive LogFC
  "#0F2080",  # unique genes with non-positive LogFC
  "#999999"   # genes with FDR above significance threshold
)

# Rename assembly identifiers to a more readable format using a named vector
assembly_map <- c(
  "GRCh37-NCBI36"       = "GRCh37/NCBI36",
  "GRCh38.12-NCBI36"    = "GRCh38.p12/NCBI36",
  "GRCh38.13-NCBI36"    = "GRCh38.p13/NCBI36",
  "CHM13v2.0-NCBI36"     = "T2T-CHM13v2.0/NCBI36",
  "GRCh38.12-GRCh37"    = "GRCh38.p12/GRCh37",
  "GRCh38.13-GRCh37"    = "GRCh38.p13/GRCh37",
  "CHM13v2.0-GRCh37"     = "T2T-CHM13v2.0/GRCh37",
  "GRCh38.13-GRCh38.12" = "GRCh38.p13/GRCh38.p12",
  "CHM13v2.0-GRCh38.12"  = "T2T-CHM13v2.0/GRCh38.p12",
  "CHM13v2.0-GRCh38.13"  = "T2T-CHM13v2.0/GRCh38.p13"
)
all_info_MSBB$assembly <- factor(
  assembly_map[all_info_MSBB$assembly],
  levels = unname(assembly_map)
)


g = ggplot(all_info_MSBB, aes(x = assembly, fill = for_plot)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(
      values = cbPalette,
      labels = c(
        "common genes LogFC > 0",
        "unique genes LogFC > 0",
        "common genes LogFC < 0",
        "unique genes LogFC < 0",
        "common genes FDR > 0.05"
      )
    ) +
    labs(
      x    = "Reference Comparison",
      y    = "Genes Expressed",
      fill = "MSBB Differential Quantification"
    ) +
    theme_minimal() +
    theme(
      axis.text.x          = element_text(angle = 45, hjust = 1, vjust = 1, size = 20),
      axis.text.y          = element_text(size = 20),
      axis.title.x.bottom  = element_text(size = 22),
      axis.title.y.left    = element_text(size = 22),
      legend.position      = c(0, 1),
      legend.justification = c(0, 1),
      legend.text          = element_text(size = 15),
      legend.title         = element_text(size = 17)
    )

ggsave(
    filename        = path_save,
    plot            = g,
    width           = 10,
    height          = 8,
    useDingbats     = FALSE
  )
