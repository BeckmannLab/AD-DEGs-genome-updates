#Figure 4A
##load in libraries
rm(list=ls()) 
library(ggplot2)
library(GGally)
library(ggplot2)
library(tidyr)
library(dplyr)

#load in ds
df_all_rosmap= readRDS("rosmap_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/rosmap_ad_interaction_results_matrix.R

#format
common = df_all_rosmap[which(!is.na(df_all_rosmap$DE_AD_sign_same)),]
df_sig = common[which(common$interaction_status!="notSignif"),]
df_sig_subset = df_sig[which(df_sig$trait=="cogdx_simplified"),]

df_sig_subset$assembly_comparison = gsub("hg19_hg18","GRCh37/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_hg19","T2T-CHM13v2.0/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_v43","T2T-CHM13v2.0/GRCh38.p13", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v30_hg19","GRCh38.p12/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_hg19","GRCh38.p13/GRCh37", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2T_hg18","T2T-CHM13v2.0/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("T2Tv30DE_cogdx_simplified","T2T-CHM13v2.0/GRCh38.p12", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v30_hg18","GRCh38.p12/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_hg18","GRCh38.p13/NCBI36", df_sig_subset$assembly_comparison)
df_sig_subset$assembly_comparison = gsub("v43_v30","GRCh38.p13/GRCh38.p12", df_sig_subset$assembly_comparison)

df_sig_subset$assembly_comparison <- factor(df_sig_subset$assembly_comparison, levels = c("GRCh37/NCBI36", "GRCh38.p12/NCBI36", "GRCh38.p13/NCBI36", "T2T-CHM13v2.0/NCBI36", "GRCh38.p12/GRCh37", "GRCh38.p13/GRCh37", "T2T-CHM13v2.0/GRCh37", "GRCh38.p13/GRCh38.p12", "T2T-CHM13v2.0/GRCh38.p12", "T2T-CHM13v2.0/GRCh38.p13"))

df_sig_subset2 = df_sig_subset[,c("assembly_comparison", "interaction_status")]

#plot
g =ggplot(df_sig_subset2, aes(x = assembly_comparison, fill = assembly_comparison)) +
  geom_bar(stat = "count", fill = "darkgreen") +  # Set fill color to dark green
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5, size = 7) +  # Add labels with count
  labs(x = "Reference pair",
       y = "# of Significant Genes") +  # Change y-axis label
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),  # Adjust size of x-axis text
        axis.text.y = element_text(size = 20),  # Adjust size of y-axis text
        axis.title = element_text(size = 20))  # Adjust size of axis titles

#save
ggsave("bar_plot_interaction_cognition.pdf",g, width = 10, height = 12,useDingbats=FALSE)
