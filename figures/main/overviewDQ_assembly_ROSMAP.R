# ==== Libraries ====
library(ggplot2)

# ==== Functions for Data Preparation & Plotting ====

# Function to set factor levels for differential-quantification categories
set_for_plot_levels <- function(df) {
  df$for_plot <- factor(
    df$for_plot,
    levels = c("New_TRUE", "New_FALSE", "Old_TRUE", "Old_FALSE", "notSignif_TRUE")
  )
  df
}

# Function to set factor levels for assembly comparisons
set_assembly_levels <- function(df) {
  df$assembly <- factor(
    df$assembly,
    levels = c(
      "GRCh37-NCBI36", "GRCh38.12-NCBI36", "GRCh38.13-NCBI36", "CHM13v2.0-NCBI36",
      "GRCh38.12-GRCh37", "GRCh38.13-GRCh37", "CHM13v2.0-GRCh37",
      "GRCh38.13-GRCh38.12", "CHM13v2.0-GRCh38.12", "CHM13v2.0-GRCh38.13"
    )
  )
  df
}

# Function to build the overview bar plot
create_overview_plot <- function(df) {
  cbPalette <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080", "#999999")

  ggplot(df, aes(x = assembly, fill = for_plot)) +
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
      fill = "ROSMAP Differential Quantification"
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
}

# Function to save a plot with dingbats disabled
save_plot_no_dingbats <- function(plot_obj, filename, width = 10, height = 8) {
  ggsave(
    filename        = filename,
    plot            = plot_obj,
    width           = width,
    height          = height,
    useDingbats     = FALSE
  )
}

# ==== Main Workflow ====

# Load data
all_info_MSBB    <- readRDS(input_rds_msbb)
info_all_ROSMAP  <- readRDS(input_rds_rosmap)
# (You can also load `all` if needed: readRDS(input_rds_all))

# Prepare data for plotting
info_all_ROSMAP <- info_all_ROSMAP %>%
  set_for_plot_levels() %>%
  set_assembly_levels()

# Build and save the plot
overview_plot <- create_overview_plot(info_all_ROSMAP)
save_plot_no_dingbats(overview_plot, output_pdf)
