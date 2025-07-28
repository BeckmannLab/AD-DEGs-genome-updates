#Supplementary Figure 7

## Clear the global environment
rm(list = ls())

library(ggplot2)
library(data.table)
############################
####DE AD results first
################################
# Clear all objects from the workspace
rm(list=ls())

# Define file paths for differential expression results for each genome assembly
path_hg18= #path to DE results
path_hg19= #path to DE results
path_v30= #path to DE results
path_v43= #path to DE results
path_T2T= #path to DE results

# Specify phenotypes or tests to process
phenos=c("braaksc_simplified","ceradsc_defvsctl","cogdx_simplified","cts_mmse30")

# List of genome assemblies
assemblies=c("hg18","hg19","v30","v43","T2T")

# Generate all unique pairs of assemblies
assembliesComb = t(combn(assemblies,2))
assembliesComb2=assembliesComb[assembliesComb[,2]=="T2T",]
assembliesComb=assembliesComb[assembliesComb[,2]!="T2T",]
assembliesComb=rbind(assembliesComb,assembliesComb2[,2:1])
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}))

# Create a grid of assembly pairs and phenotypes
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
grid$Var4[grid$Var4=="ceradsc_defvsctl"]="ceradsc_test.txt"
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")

# Prepare paths for each assembly
grid3=as.data.frame(assemblies)
grid3$assemblies_path=paste("path",grid3$assemblies, sep = "_")

# Loop through each assembly to read and merge DE files
for(i in 1:nrow(grid3)){
path=get(grid3[i,"assemblies_path"])
DE_files=list.files(path,pattern="^DE")
all=data.table()
    for(DE_file in DE_files){
        tmp=fread(paste0(path,DE_file))
        phenotype=unlist(lapply(strsplit(DE_file,"_"),function(x){paste0(x[20],"_",x[21])}))
        colnames(tmp)[2:ncol(tmp)]=paste0(colnames(tmp)[2:ncol(tmp)],"|",phenotype)
        if(nrow(all)<1){
            all=tmp
        }else{
            all=merge(all,tmp,by="V1",all=TRUE)
        }}
colnames(all)=paste0(colnames(all),"_",grid3[i,"assemblies"])
colnames(all)[1]="gene_id"
value_assigned=as.character(grid3[i, "assemblies"])
assign(value_assigned,all)
}

# Merge results from all assemblies into one data frame
all_results_AD=merge(hg18,hg19,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,v30,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,v43,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,T2T,by="gene_id",all=TRUE)
all_results_AD=as.data.frame(all_results_AD)

# Reinitialize combinations for significance categorization
assembliesComb = t(combn(assemblies,2))
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[2],x[1])}))

# Recreate grid for significance loop
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
grid$Var4[grid$Var4=="ceradsc_defvsctl"]="ceradsc_test.txt"
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")

# Copy merged results for analysis
all_results = all_results_AD

# Loop to categorize significance across assembly pairs for each phenotype
for(i in 1:nrow(grid)){
    # Define output column for this comparison
    outCol=paste0(grid[i,"assembly2"],grid[i,"assembly1"],"DE_",grid[i,"pheno"])
    # Identify adjusted p-value columns for each assembly
    colOfInterest1=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly1"])
    colOfInterest2=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly2"])
    # Initialize output column
    all_results[,outCol]=""

    # Default label when both present but not significant
    all_results[!is.na(all_results[,colOfInterest1]) & !is.na(all_results[,colOfInterest2]),outCol]="notSignif"

    # Cases where one is NA and the other significant
    all_results[which(is.na(all_results[, colOfInterest1]) & all_results[, colOfInterest2] <= 0.05), outCol] = paste0(grid[i,"assembly2"],"Signif")
    all_results[which(is.na(all_results[,colOfInterest2]) & all_results[,colOfInterest1]<=0.05),outCol]=paste0(grid[i,"assembly1"],"Signif")

    # Both significant
    all_results[which(all_results[,colOfInterest1]<=0.05 & all_results[,colOfInterest2]<=0.05),outCol]="bothSignif"

    # One significant, one not
    all_results[which(all_results[,colOfInterest1]<=0.05 & all_results[,colOfInterest2]>0.05),outCol]=paste0(grid[i,"assembly1"],"Signif")
    all_results[which(all_results[,colOfInterest1]>0.05 & all_results[,colOfInterest2]<=0.05),outCol]=paste0(grid[i,"assembly2"],"Signif")

    # Cases where one is NA but not significant
    all_results[which(is.na(all_results[,colOfInterest1]) & all_results[,colOfInterest2]>0.05),outCol]="notSignif"
    all_results[which(is.na(all_results[,colOfInterest2]) & all_results[,colOfInterest1]>0.05),outCol]="notSignif"

    # Replace empty strings with NA
    all_results[all_results[,outCol]=="",outCol]=NA

    # Set factor levels for categorization
    all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly1"],"Signif"),paste0(grid[i,"assembly2"],"Signif"),"bothSignif","notSignif"))
}

# Convert to data frame and save results
all_results=as.data.frame(all_results)

# 2) Define comparisons and assembly IDs
comparisons <- c(
  "NCBI36 vs GRCh37",      "NCBI36 vs GRCh38.p12", "NCBI36 vs GRCh38.p13",
  "NCBI36 vs CHM13v2.0","GRCh37 vs GRCh38.p12",  "GRCh37 vs GRCh38.p13",
  "GRCh37 vs CHM13v2.0","GRCh38.p12 vs GRCh38.p13",
  "GRCh38.p12 vs CHM13v2.0","GRCh38.p13 vs CHM13v2.0"
)
older_ids <- c("hg18","hg18","hg18","hg18","hg19","hg19","hg19","v30","v30","v43")
newer_ids <- c("hg19","v30","v43","T2T","v30","v43","T2T","v43","T2T","T2T")

# 3) Build the plotting data.frame
plot_data <- do.call(rbind, Map(function(comp, old, new) {
  xcol  <- paste0("logFC|ceradsc_test.txt_", old)
  ycol  <- paste0("logFC|ceradsc_test.txt_", new)
  decol <- paste0(new, old, "DE_ceradsc_defvsctl")

  status <- ifelse(
    all_results[[decol]] == "notSignif",                "Neither reference significant",
  ifelse(all_results[[decol]] == paste0(old, "Signif"),  "Older reference significant",
  ifelse(all_results[[decol]] == paste0(new, "Signif"),  "Newer reference significant",
  ifelse(all_results[[decol]] == "bothSignif",           "Both references significant",
                                                            NA))))

  data.frame(
    Comparison = comp,
    x          = all_results[[xcol]],
    y          = all_results[[ycol]],
    DE_status  = factor(
      status,
      levels = c(
        "Neither reference significant",
        "Both references significant",
        "Older reference significant",
        "Newer reference significant"
      )
    ),
    stringsAsFactors = FALSE
  )
}, comparisons, older_ids, newer_ids))

# 4) Fix facet order
desired_order <- c(
  "NCBI36 vs GRCh37",
  "NCBI36 vs GRCh38.p12",
  "NCBI36 vs GRCh38.p13",
  "NCBI36 vs CHM13v2.0",
  "GRCh37 vs GRCh38.p12",
  "GRCh37 vs GRCh38.p13",
  "GRCh37 vs CHM13v2.0",
  "GRCh38.p12 vs GRCh38.p13",
  "GRCh38.p12 vs CHM13v2.0",
  "GRCh38.p13 vs CHM13v2.0"
)
plot_data$Comparison <- factor(plot_data$Comparison, levels = desired_order)
plot_data$DE_status <- factor(plot_data$DE_status, levels = c("Newer reference significant", "Older reference significant", "Both references significant", "Neither reference significant"))

# 5) Plot: first draw Neither & Both, then Older & Newer on top
g <- ggplot(plot_data, aes(x, y)) +
  # layer 1: neither (grey) + both (green)
  geom_point(
    data = subset(plot_data, DE_status %in% c("Neither reference significant", "Both references significant")),
    aes(color = DE_status),
    size  = 2.5,
    alpha = 0.4
  ) +
  # layer 2: older (red) + newer (blue) ON TOP
  geom_point(
    data = subset(plot_data, DE_status %in% c("Older reference significant", "Newer reference significant")),
    aes(color = DE_status),
    size  = 2.5,
    alpha = 0.4
  ) +
  facet_wrap(~Comparison, nrow = 2) +
  # scale_x_continuous(expand = expansion(add = 0.05)) +
  # scale_y_continuous(expand = expansion(add = 0.05)) +
  scale_color_manual(
    name   = "AD DE Status",
    # specify breaks here to force legend order
    breaks = c(
      "Older reference significant",
      "Newer reference significant",
      "Both references significant",
      "Neither reference significant"
    ),
    values = c(
      "Older reference significant"    = "red",
      "Newer reference significant"    = "blue",
      "Neither reference significant"  = "grey80",
      "Both references significant"    = "darkgreen"
    ),
    guide = guide_legend(
      override.aes   = list(size = 5),
      direction      = "horizontal",
      title.position = "left",
      title.hjust    = 0
    )
  ) +
  labs(
    title = "CERAD",
    x     = "DE LogFC older reference",
    y     = "DE LogFC newer reference"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title         = element_text(size = 20, face = "bold", hjust = 0.5),
    strip.text         = element_text(size = 14),
    axis.title         = element_text(size = 16),
    axis.text          = element_text(size = 14),
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.box         = "horizontal",
    legend.title       = element_text(size = 20),
    legend.text        = element_text(size = 18),
    legend.title.align = 0,
    legend.key.size    = unit(1.5, "lines"),
    legend.background  = element_rect(color = "black", fill = "white"),
    panel.spacing      = unit(0.5, "lines"),
    plot.margin        = margin(10, 10, 10, 10)
  ) +
  coord_fixed()  # placed last so no subsequent scale/facet call can override it


ggsave("cerad_scatter_plot.pdf",g, useDingbats = FALSE, width = 17.5, height= 9)
