#Supplementary Figure 12
# ---------------------------------
# Load Input Data Tables
# ---------------------------------
combo_var <- read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/combo_var.txt", header = F, sep = "\t")
# structure(list(V1 = 1:10, V2 = c("hg19_hg18", "v30_hg18", "v43_hg18", 
# "T2T_hg18", "v30_hg19", "v43_hg19", "T2T_hg19", "v30_v43", "T2T_v30", 
# "T2T_v43")), class = "data.frame", row.names = c(NA, -10L))

all_results = readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/analysis/msbb_all_resultsDE_5.02.24.RDS") #see AD-DEGs-genome-updates/misc/msbb_all_resultsDE_5.02.24.R

var = read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/analysis/var.txt", 
  header = F, sep = "\t")
# structure(list(V1 = c("hg18hg19DE_ceradsc_defvsctl", "hg18v30DE_ceradsc_defvsctl", 
# "hg18v43DE_ceradsc_defvsctl", "hg19v30DE_ceradsc_defvsctl", "hg19v43DE_ceradsc_defvsctl", 
# "v30v43DE_ceradsc_defvsctl", "T2Thg18DE_ceradsc_defvsctl", "T2Thg19DE_ceradsc_defvsctl", 
# "T2Tv30DE_ceradsc_defvsctl", "T2Tv43DE_ceradsc_defvsctl"), V2 = c("hg18Signif", 
# "hg18Signif", "hg18Signif", "hg19Signif", "hg19Signif", "v30Signif", 
# "hg18Signif", "hg19Signif", "v30Signif", "v43Signif"), V3 = c("hg19Signif", 
# "v30Signif", "v43Signif", "v30Signif", "v43Signif", "v43Signif", 
# "T2TSignif", "T2TSignif", "T2TSignif", "T2TSignif")), class = "data.frame", row.names = c(NA, 
# -10L))

# -----------------------------
# Load Required R Libraries
# -----------------------------
library(tidyr)    
library(ggplot2)  
library(GGally)  
library(patchwork) 

# -----------------------------------
# Rename Assembly Comparisons
# -----------------------------------
combo_var[1,1]  = "GRCh37/NCBI36"
combo_var[2,1]  = "GRCh38.p12/NCBI36"
combo_var[3,1]  = "GRCh38.p13/NCBI36"
combo_var[4,1]  = "CHM13v2.0/NCBI36"
combo_var[5,1]  = "GRCh38.p12/GRCh37"
combo_var[6,1]  = "GRCh38.p13/GRCh37"
combo_var[7,1]  = "CHM13v2.0/GRCh37"
combo_var[8,1]  = "GRCh38.p13/GRCh38.p12"
combo_var[9,1]  = "CHM13v2.0/GRCh38.p12"
combo_var[10,1] = "CHM13v2.0/GRCh38.p13"

# --------------------------------
# Prepare Data for Plotting
# --------------------------------

CDR_simplified = all_results[, c(
  "hg18DE_CDR_simplified", 
  "hg19DE_CDR_simplified", 
  "v30DE_CDR_simplified", 
  "v43DE_CDR_simplified", 
  "T2TDE_CDR_simplified"
)]

CERJ_defvsctl = all_results[, c(
  "hg18DE_CERJ_defvsctl", 
  "hg19DE_CERJ_defvsctl", 
  "v30DE_CERJ_defvsctl", 
  "v43DE_CERJ_defvsctl", 
  "T2TDE_CERJ_defvsctl"
)]

PlaqueMean_test.txt = all_results[, c(
  "hg18DE_PlaqueMean_test.txt", 
  "hg19DE_PlaqueMean_test.txt", 
  "v30DE_PlaqueMean_test.txt", 
  "v43DE_PlaqueMean_test.txt", 
  "T2TDE_PlaqueMean_test.txt"
)]

CERJ_defvsctl_results <- gather(
  CERJ_defvsctl, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE
)
CERJ_defvsctl_results$test = "CERJ_defvsctl"

CDR_simplified_results <- gather(
  CDR_simplified, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE
)
CDR_simplified_results$test = "CDR_simplified"

PlaqueMean_test.txt_results <- gather(
  PlaqueMean_test.txt, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE
)
PlaqueMean_test.txt_results$test = "PlaqueMean_test.txt"

data_long= rbind(CERJ_defvsctl_results,CDR_simplified_results,PlaqueMean_test.txt_results)

# -----------------------
# Build the Plot
# -----------------------
a = ggplot(data_long, aes(x = Category, fill = Count)) +
  geom_bar(stat = "count", position = "stack") +
  xlab("Reference") +
  ylab("Genes Expressed") +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 20, hjust = 1, angle = 90),
    axis.text.y  = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text  = element_text(size = 17),
    legend.title = element_text(size = 18),
    plot.title   = element_text(size = 20, face = "bold", hjust = 0.5)
  ) +
  scale_fill_manual(
    values = c("grey", "darkgreen"),
    labels = expression("" > 0.05, "" <= 0.05)
  ) +
  labs(fill = "FDR") +
  geom_text(
    stat    = "count",
    aes(label = ..count..),
    vjust   = -0.5,
    size    = 3
  ) +
  scale_x_discrete(labels = c("NCBI36", "GRCh37", "GRCh38.p12", "GRCh38.p13", "CHM13v2.0")) +
  facet_wrap(
    ~ test,
    scales = "free_x",
    labeller = as_labeller(c(
      CERJ_defvsctl           = "CERJ",
      CDR_simplified     = "CDR",
      PlaqueMean_test.txt   = "Plaque Mean"
    ))
  )


# --------------------------------
# Save Plot to Files
# --------------------------------

ggsave(a, filename    = paste0("msbb_overview_DE_all_tests.pdf"),useDingbats = FALSE, width = 10, height = 6)

