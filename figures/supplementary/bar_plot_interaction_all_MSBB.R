#Supplementary figure 20

##load in libraries
rm(list=ls()) 
library(ggplot2)
library(GGally)
library(ggplot2)
library(tidyr)
library(dplyr)

#load in ds
df_all_msbb = readRDS("MSBB_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/MSBB_ad_interaction_results_matrix.R

#format
common = df_all_msbb[which(!is.na(df_all_msbb$DE_AD_sign_same)),]
df_sig = common[which(common$interaction_status!="notSignif"),]
df_sig_subset <- df_sig[df_sig$trait %in% c("CDR_simplified", "PlaqueMean", "CERJ_defvsctl"), ]

df_sig_subset$assembly_comparison = gsub("hg19_hg18","GRCh37/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_hg19","CHM13v2.0/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_v43","CHM13v2.0/GRCh38.p13", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v30_hg19","GRCh38.p12/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_hg19","GRCh38.p13/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_hg18","CHM13v2.0/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_v30","CHM13v2.0/GRCh38.p12", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v30_hg18","GRCh38.p12/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_hg18","GRCh38.p13/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_v30","GRCh38.p13/GRCh38.p12", df_sig_subset$assembly_comparison)

df_sig_subset$assembly_comparison <- factor(df_sig_subset$assembly_comparison, levels = c("GRCh37/NCBI36", "GRCh38.p12/NCBI36", "GRCh38.p13/NCBI36", "CHM13v2.0/NCBI36", "GRCh38.p12/GRCh37", "GRCh38.p13/GRCh37", "CHM13v2.0/GRCh37", "GRCh38.p13/GRCh38.p12", "CHM13v2.0/GRCh38.p12", "CHM13v2.0/GRCh38.p13"))

df_sig_subset2 = df_sig_subset[,c("assembly_comparison", "interaction_status", "trait")]
df_sig_subset2$new_name = gsub("CDR_simplified", "CDR", df_sig_subset2$trait)
df_sig_subset2$new_name = gsub("PlaqueMean", "Plaque mean density", df_sig_subset2$trait)

#plot
g <- ggplot(df_sig_subset2, aes(x = assembly_comparison, fill = assembly_comparison)) +
  geom_bar(stat = "count", fill = "darkgreen") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 7) +
  labs(
    x = "Reference pair",
    y = "# of Significant Genes"
  ) +
  facet_wrap(~ new_name, scales = "free_y", ncol = 4) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +  # remove excess space
  coord_cartesian(ylim = c(0, NA)) +  # ensures no padding above max bar
  theme_minimal() +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1, size = 20),
    axis.text.y   = element_text(size = 20),
    axis.title    = element_text(size = 20),
    panel.spacing = unit(0.3, "lines")
  )

#save
ggsave("bar_plot_interaction_all_MSBB.pdf",g, width = 22, height = 11,useDingbats=FALSE)
