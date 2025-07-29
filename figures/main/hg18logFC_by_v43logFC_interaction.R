#Figure 5B

library(ggrepel)
library(ggplot2)

interaction_rosmap=readRDS("interaction_andAD_DE_with_symbols_4.10.24.RDS") #see AD-DEGs-genome-updates/misc/interaction_andAD_DE_with_symbols_4.10.24.RDS
mapping=readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.RDS

var=read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/analysis/var.txt", header = F, sep="\t")
# structure(list(V1 = c("hg18hg19DE_ceradsc_defvsctl", "hg18v30DE_ceradsc_defvsctl", 
# "hg18v43DE_ceradsc_defvsctl", "hg19v30DE_ceradsc_defvsctl", "hg19v43DE_ceradsc_defvsctl", 
# "v30v43DE_ceradsc_defvsctl", "T2Thg18DE_ceradsc_defvsctl", "T2Thg19DE_ceradsc_defvsctl", 
# "T2Tv30DE_ceradsc_defvsctl", "T2Tv43DE_ceradsc_defvsctl"), V2 = c("hg18Signif", 
# "hg18Signif", "hg18Signif", "hg19Signif", "hg19Signif", "v30Signif", 
# "hg18Signif", "hg19Signif", "v30Signif", "v43Signif"), V3 = c("hg19Signif", 
# "v30Signif", "v43Signif", "v30Signif", "v43Signif", "v43Signif", 
# "T2TSignif", "T2TSignif", "T2TSignif", "T2TSignif")), class = "data.frame", row.names = c(NA, 
# -10L))

#formatting
colnames(interaction_rosmap) <- gsub("\\|", "_", colnames(interaction_rosmap))
interaction_rosmap=as.data.frame(interaction_rosmap)
interaction_rosmap$for_plot_v43hg18=ifelse(interaction_rosmap$interaction_v43hg18DE_cogdx_simplified == "notSignif", "notSignif", "Signif")
interaction_rosmap$for_plot_v43hg18 <- factor(interaction_rosmap$for_plot_v43hg18, levels = c("Signif", "notSignif"))
interaction_rosmap=na.omit(interaction_rosmap)

#plot
g = ggplot(interaction_rosmap[order(interaction_rosmap$for_plot_v43hg18, decreasing = TRUE),],
           aes(x = DEAD_logFC_cogdx_simplified_v43, 
               y = DEAD_logFC_cogdx_simplified_hg18,
               color = for_plot_v43hg18, 
               shape = DEAD_v43hg18DE_cogdx_simplified)) + 
  geom_point(alpha = 0.5, size = 4) + 
  geom_text_repel(data = subset(interaction_rosmap, gene_symbol %in% c("RNF130")),
                  mapping = aes(label = gene_symbol),
                  min.segment.length = 0,
                  box.padding = 0.5,
                  nudge_y = 0.1, 
                  show.legend = FALSE, 
                  size = 10) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c(notSignif = "gray", Signif = "darkgreen"),
                     labels = expression(paste("" <= 0.05), paste("" > 0.05))) + 
  scale_shape_manual(values = c(16, 17, 15, 3),
                     labels = expression(paste("" <= 0.05, '    for NCBCI36'),
                                          paste("" <= 0.05, '   for GRCh38.p13'),
                                          paste("" <= 0.05, '   for both'), 
                                          paste("" > 0.05))) + 
  coord_fixed() + 
  xlab("GRCh38.p13 LogFC") + 
  ylab("NCBI36 LogFC") + 
  labs(color = 'FDR for DE AD', shape = "FDR for Interaction") +
  guides(color = guide_legend(override.aes = list(shape = 4))) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title = element_text(size = 28),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 23),
    # legend.background = element_rect(color = "black", fill = "white", linewidth = 0.8),
    # legend.box.margin = margin(6, 6, 6, 6)
  )

ggsave("hg18logFC_by_v43logFC_interaction.pdf",g, width = 13, height = 10,useDingbats=FALSE)


