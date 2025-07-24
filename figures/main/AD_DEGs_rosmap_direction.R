#Figure 3D

rm(list=ls())
library(ggplot2)
library(data.table)
library(corrplot)
library(reshape2)
library(dplyr)
library(tidyr)

##for loop for plotting 
combo_var <- read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/combo_var.txt", header = F, sep="\t")
# structure(list(V1 = 1:10, V2 = c("hg19_hg18", "v30_hg18", "v43_hg18", 
# "T2T_hg18", "v30_hg19", "v43_hg19", "T2T_hg19", "v30_v43", "T2T_v30", 
# "T2T_v43")), class = "data.frame", row.names = c(NA, -10L))


overview <- read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/old_Rosmap_combos/overview_assemblyDE.txt", sep="\t")
# structure(list(assembly = c("NCBI36_GRCh37", "NCBI36_GRCh38.12", 
# "NCBI36_GRCh38.13", "NCBI36_CHM13v2.0", "GRCh37_GRCh38.12", "GRCh37_GRCh38.13", 
# "GRCh37_CHM13v2.0", "GRCh38.12_GRCh38.13", "GRCh38.12_CHM13v2.0", 
# "GRCh38.13_CHM13v2.0"), genes = c(14627L, 14306L, 14284L, 14265L, 
# 17476L, 17352L, 17082L, 19220L, 18595L, 18916L), non_DE = c(95L, 
# 127L, 85L, 83L, 209L, 204L, 200L, 216L, 269L, 264L), DE = c(14532L, 
# 14179L, 14199L, 14182L, 17267L, 17148L, 16882L, 19004L, 18326L, 
# 18652L)), class = "data.frame", row.names = c(NA, -10L))

combo_var[1,1]="GRCh37/NCBI36"
combo_var[2,1]="GRCh38.p12/NCBI36"
combo_var[3,1]="GRCh38.p13/NCBI36"
combo_var[4,1]="CHM13v2.0/NCBI36"
combo_var[5,1]="GRCh38.p12/GRCh37"
combo_var[6,1]="GRCh38.p13/GRCh37"
combo_var[7,1]="CHM13v2.0/GRCh37"
combo_var[8,1]="GRCh38.p13/GRCh38.p12"
combo_var[9,1]="CHM13v2.0/GRCh38.p12"
combo_var[10,1]="CHM13v2.0/GRCh38.p13"

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
var = as.data.frame(var)

var$V1=gsub("_ceradsc_defvsctl", "_cogdx_simplified",var$V1)
var$V1= gsub("hg18hg19DE", "hg19hg18DE",var$V1)
var$V1= gsub("hg18v30DE", "v30hg18DE",var$V1)
var$V1= gsub("hg18v43DE", "v43hg18DE",var$V1)
var$V1= gsub("hg19v30DE", "v30hg19DE",var$V1)
var$V1= gsub("hg19v43DE", "v43hg19DE",var$V1)
var$V1= gsub("v30v43DE", "v43v30DE",var$V1)

############################
####DE AD results first
################################

# Define file paths for differential expression results for each genome assembly
path_hg18= #path to DE results
path_hg19= #path to DE results
path_v30= #path to DE results
path_v43= #path to DE results
path_T2T= #path to DE results

path_hg18="/sc/arion/projects/mscic1/results/anina/fun_project_4.23/RosmapHG18/DE/"
path_hg19="/sc/arion/projects/mscic1/results/anina/fun_project_4.23/RosmapHG19/DE/"
path_v30="/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmapv30/DE/"
path_v43="/sc/arion/projects/mscic1/results/anina/fun_project_4.23/Rosmapv43/DE/"
path_T2T="/sc/arion/projects/mscic1/results/anina/fun_project_4.23/RosmapT2T/DE/"

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


cogdx_subset=all_results[,grep("cogdx_simplified",colnames(all_results))]

overview[,1]=combo_var[,1]
overview <- as.data.frame(overview[,c(1)])
for (i in seq_len(nrow(var))) {
  x <- cogdx_subset[[var[i,1]]]
  older=length(which(x==var[i,2]))
  newer=length(which(x==var[i,3]))
  both=length(which(x=="bothSignif"))
  overview[i,2] <- older
  overview[i,3] <- newer
  overview[i,4] <- both
}

colnames(overview)=c("assembly", "older", "newer", "both")

overview_t <- pivot_longer(overview, c(older, newer, both))


overview_t$assembly <- factor(overview_t$assembly, levels = c("GRCh37/NCBI36", "GRCh38.p12/NCBI36", "GRCh38.p13/NCBI36", "CHM13v2.0/NCBI36", "GRCh38.p12/GRCh37", "GRCh38.p13/GRCh37", "CHM13v2.0/GRCh37", "GRCh38.p13/GRCh38.p12", "CHM13v2.0/GRCh38.p12", "CHM13v2.0/GRCh38.p13"))

overview_t$name = gsub("older", "Older reference", overview_t$name)
overview_t$name = gsub("newer", "Newer reference", overview_t$name)
overview_t$name = gsub("both", "Both references", overview_t$name)

overview_t$name = factor(overview_t$name, levels = c("Older reference", "Newer reference", "Both references"))

g <- ggplot(overview_t, aes(fill = name, y = value, x = assembly)) + 
  geom_bar(position = "dodge", stat = "identity", width = 0.7) + 
  scale_fill_manual(
    breaks = c("Older reference", "Newer reference", "Both references"),
    values = c(
      "Older reference"  = "darkblue",
      "Newer reference"  = "lightblue",
      "Both references"  = "darkgreen"
    ),
    labels = c("Older reference", "Newer reference", "Both references")
  ) +
  labs(x = "Reference Comparison", y = "Number of DEGs") +
  guides(fill = guide_legend(title = "AD DEG in")) +
  theme_bw(base_size = 14) +
  theme(
    axis.text          = element_text(size = 22),
    axis.title         = element_text(size = 24),
    plot.title         = element_text(hjust = 0.5, size = 16),
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 19),
    legend.title       = element_text(size = 20),
    legend.text        = element_text(size = 18),
    legend.background  = element_rect(fill = "white", colour = "black", size = 0.5),
    legend.key         = element_rect(fill = NA, colour = NA),
    legend.position    = c(0.02, 0.98),
    legend.justification = c(0, 1)
  )

ggsave(g, filename=paste0("AD_DEGs_rosmap_direction.pdf"), width = 12, height = 9, useDingbats=FALSE)


