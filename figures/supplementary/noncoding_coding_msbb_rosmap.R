#Supplementary figure 5C
#--------------------------------------------
# Files
#--------------------------------------------
rosmap_data_path <- "with_coding_status_all_info_ROSMAP.RDS" #see AD-DEGs-genome-updates/misc/with_coding_status_all_info_ROSMAP.R
msbb_data_path   <- "with_coding_status_all_info_MSBB.RDS" #see AD-DEGs-genome-updates/misc/with_coding_status_all_info_MSBB.R
#--------------------------------------------
# Load libraries
#--------------------------------------------
library(ggplot2)
library(data.table)
library(tidyr)
library(GGally)
library(stringr)
library(patchwork)

#--------------------------------------------
# Constants and mappings
#--------------------------------------------
assembly_mapping <- c(
  "GRCh37-NCBI36"       = "GRCh37/NCBI36",
  "GRCh38.12-NCBI36"    = "GRCh38.p12/NCBI36",
  "GRCh38.13-NCBI36"    = "GRCh38.p13/NCBI36",
  "CHM13v2.0-NCBI36"    = "T2T-CHM13v2.0/NCBI36",
  "GRCh38.12-GRCh37"    = "GRCh38.p12/GRCh37",
  "GRCh38.13-GRCh37"    = "GRCh38.p13/GRCh37",
  "CHM13v2.0-GRCh37"    = "T2T-CHM13v2.0/GRCh37",
  "GRCh38.13-GRCh38.12" = "GRCh38.p13/GRCh38.p12",
  "CHM13v2.0-GRCh38.12" = "T2T-CHM13v2.0/GRCh38.p12",
  "CHM13v2.0-GRCh38.13" = "T2T-CHM13v2.0/GRCh38.p13"
)
assembly_levels  <- unname(assembly_mapping)
levels_for_plot  <- c("New_TRUE", "New_FALSE", "Old_TRUE", "Old_FALSE", "notSignif_TRUE")
cbPalette        <- c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080", "#999999")

variable_names <- c(
  protein_coding = "Protein Coding",
  non_coding     = "Non Coding",
  New            = "Quantified higher\nin newer assembly",
  Old            = "Quantified higher\nin older assembly"
)
variable_labeller <- function(variable, value) variable_names[value]

#--------------------------------------------
# Helper functions
#--------------------------------------------
clean_assembly_names <- function(df) {
  df$assembly <- str_replace_all(df$assembly, assembly_mapping)
  df$assembly <- factor(df$assembly, levels = assembly_levels)
  df$for_plot  <- factor(df$for_plot, levels = levels_for_plot)
  return(df)
}

generate_plot <- function(data, title_text) {
  ggplot(data, aes(x = assembly,
                   y = after_stat(prop),
                   fill = for_plot,
                   by = assembly)) +
    geom_bar(position = "stack", stat = "prop") +
    theme_minimal(base_size = 16) +
    labs(
      x     = "Reference",
      y     = "Proportion of Genes Expressed",
      title = title_text
    ) +
    theme(
      axis.text.x     = element_text(angle = 90, size = 16),
      axis.text.y     = element_text(size = 16),
      axis.title.x    = element_text(size = 18),
      axis.title.y    = element_text(size = 18),
      plot.title      = element_text(hjust = 0.5, size = 24),
      legend.text     = element_text(size = 18),
      legend.title    = element_blank(),
      strip.text      = element_text(size = 18),
      legend.position = "bottom"
    ) +
    scale_fill_manual(
      values = cbPalette,
      labels = c(
        'common genes LogFC > 0',
        'unique genes LogFC > 0',
        'common genes LogFC < 0',
        'unique genes LogFC < 0',
        'common genes FDR > 0.05'
      )
    ) +
    facet_grid(status ~ direction, labeller = variable_labeller) +
    geom_hline(yintercept = 0.5, color = "white", size = 0.5)
}

#--------------------------------------------
# Main script execution
#--------------------------------------------
rosmap <- readRDS(rosmap_data_path)
rosmap <- clean_assembly_names(rosmap)
g     <- generate_plot(rosmap, "ROSMAP")

msbb  <- readRDS(msbb_data_path)
msbb  <- clean_assembly_names(msbb)
b     <- generate_plot(msbb, "MSBB")

w <- g + b + plot_layout(guides = "collect") &
     theme(legend.position = 'bottom')

ggsave(
  filename    = file.path(output_dir, "noncoding_coding_msbb_rosmap.pdf"),
  plot        = w,
  useDingbats = FALSE
)
