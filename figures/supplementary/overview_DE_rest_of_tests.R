#Supplementary Figure 6
# ---------------------------------
# Load Input Data Tables
# ---------------------------------
combo_var <- read.delim(
  "/sc/arion/projects/mscic1/results/anina/fun_project_4.23/combo_var.txt", header = F, sep = "\t")

# structure(list(V1 = c("GRCh37/NCBI36", "GRCh38.p12/NCBI36", "GRCh38.p13/NCBI36", 
# "CHM13v2.0/NCBI36", "GRCh38.p12/GRCh37", "GRCh38.p13/GRCh37", 
# "CHM13v2.0/GRCh37", "GRCh38.p13/GRCh38.p12", "CHM13v2.0/GRCh38.p12", 
# "CHM13v2.0/GRCh38.p13"), V2 = c("hg19_hg18", "v30_hg18", "v43_hg18", 
# "T2T_hg18", "v30_hg19", "v43_hg19", "T2T_hg19", "v30_v43", "T2T_v30", 
# "T2T_v43")), row.names = c(NA, -10L), class = "data.frame")

overview <- read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/old_Rosmap_combos/overview_assemblyDE.txt",sep = "\t")

# structure(list(assembly = c("NCBI36_GRCh37", "NCBI36_GRCh38.12", 
# "NCBI36_GRCh38.13", "NCBI36_CHM13v2.0", "GRCh37_GRCh38.12", "GRCh37_GRCh38.13", 
# "GRCh37_CHM13v2.0", "GRCh38.12_GRCh38.13", "GRCh38.12_CHM13v2.0", 
# "GRCh38.13_CHM13v2.0"), genes = c(14627L, 14306L, 14284L, 14265L, 
# 17476L, 17352L, 17082L, 19220L, 18595L, 18916L), non_DE = c(95L, 
# 127L, 85L, 83L, 209L, 204L, 200L, 216L, 269L, 264L), DE = c(14532L, 
# 14179L, 14199L, 14182L, 17267L, 17148L, 16882L, 19004L, 18326L, 
# 18652L)), class = "data.frame", row.names = c(NA, -10L))

all_results = readRDS("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/analysis/all_resultsDE_7.24.25.RDS") #see AD-DEGs-genome-updates/misc/all_resultsDE_7.24.25.R

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
glibrary(GGally)  
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
braaksc_simplified = all_results[, c(
  "hg18DE_braaksc_simplified", 
  "hg19DE_braaksc_simplified", 
  "v30DE_braaksc_simplified", 
  "v43DE_braaksc_simplified", 
  "T2TDE_braaksc_simplified"
)]

ceradsc_defvsctl = all_results[, c(
  "hg18DE_ceradsc_defvsctl", 
  "hg19DE_ceradsc_defvsctl", 
  "v30DE_ceradsc_defvsctl", 
  "v43DE_ceradsc_defvsctl", 
  "T2TDE_ceradsc_defvsctl"
)]

cts_mmse30 = all_results[, c(
  "hg18DE_cts_mmse30", 
  "hg19DE_cts_mmse30", 
  "v30DE_cts_mmse30", 
  "v43DE_cts_mmse30", 
  "T2TDE_cts_mmse30"
)]

ceradsc_defvsctl_results <- gather(
  ceradsc_defvsctl, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE  
)
ceradsc_defvsctl_results$test = "ceradsc_defvsctl"

braaksc_simplified_results <- gather(
  braaksc_simplified, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE  
)
braaksc_simplified_results$test = "braaksc_simplified"

cts_mmse30_results <- gather(
  cts_mmse30, 
  key   = "Category", 
  value = "Count",
  na.rm = TRUE  
)
cts_mmse30_results$test = "cts_mmse30"

data_long= rbind(ceradsc_defvsctl_results, braaksc_simplified_results,cts_mmse30_results)

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
      cts_mmse30           = "MMSE30",
      ceradsc_defvsctl     = "CERAD",
      braaksc_simplified   = "Braak"
    ))
  )

# --------------------------------
# Save Plot to Files
# --------------------------------

ggsave(a, filename    = paste0(
    "overview_DE_rest_of_tests.pdf"),useDingbats = FALSE, width = 10, height = 6)
