#Supplementary figure 3C

#libraries
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# function
make_logFC_plot <- function(data, title) {
  ggplot(data, aes(x = logFC)) +
    geom_density(fill = "skyblue", color = "darkblue", alpha = 0.7) +
    labs(title = title, x = "LogFC", y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
}

#load in countrs
hg19_hg18=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v30_hg18=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v30_hg19=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v30_v43=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")  #counts
v43_hg18=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v43_hg19=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
T2T_hg18=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
T2T_hg19=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
T2T_v30=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
T2T_v43=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts

#formatting
hg19_hg18$assembly = "GRCh37-NCBI36"
v30_hg18$assembly = "GRCh38.p12-NCBI36"
v30_hg19$assembly = "GRCh38.p12-GRCh37"
v30_v43$assembly = "GRCh38.p13-GRCh38.p12"
v43_hg18$assembly = "GRCh38.p13-NCBI36"
v43_hg19$assembly = "GRCh38.p13-GRCh37"
T2T_hg18$assembly = "CHM13v2.0-NCBI36"
T2T_hg19$assembly = "CHM13v2.0-GRCh37"
T2T_v30$assembly = "CHM13v2.0-GRCh38.p12"
T2T_v43$assembly = "CHM13v2.0-GRCh38.p13"

#combine
all= rbind(hg19_hg18, v30_hg18,v30_hg19, v30_v43,v43_hg18,v43_hg19,T2T_hg18, T2T_hg19,T2T_v30, T2T_v43)

all$assembly <- factor(all$assembly , levels = c("GRCh37-NCBI36","GRCh38.p12-NCBI36", "GRCh38.p13-NCBI36", "CHM13v2.0-NCBI36", "GRCh38.p12-GRCh37", "GRCh38.p13-GRCh37", "CHM13v2.0-GRCh37", "GRCh38.p13-GRCh38.p12", "CHM13v2.0-GRCh38.p12", "CHM13v2.0-GRCh38.p13"))

#plot
my_colors <- brewer.pal(10, "Paired")

# List of data-title pairs
plot_list <- list(
  make_logFC_plot(hg19_hg18, "LogFC GRCh37-NCBI36"),
  make_logFC_plot(v30_hg18, "LogFC GRCh38.p12-NCBI36"),
  make_logFC_plot(v43_hg18, "LogFC GRCh38.p13-NCBI36"),
  make_logFC_plot(T2T_hg18, "LogFC CHM13v2.0-NCBI36"),
  make_logFC_plot(v30_hg19, "LogFC GRCh38.p12-GRCh37"),
  make_logFC_plot(v43_hg19, "LogFC GRCh38.p13-GRCh37"),
  make_logFC_plot(T2T_hg19, "LogFC CHM13v2.0-GRCh37"),
  make_logFC_plot(v30_v43, "LogFC GRCh38.p13-GRCh38.p12"),
  make_logFC_plot(T2T_v30, "LogFC CHM13v2.0-GRCh38.p12"),
  make_logFC_plot(T2T_v43, "LogFC CHM13v2.0-GRCh38.p13")
)

# Export to PDF
pdf("logFC_distribution.pdf")
grid.arrange(grobs = plot_list, ncol = 2, nrow = 5)
dev.off()
