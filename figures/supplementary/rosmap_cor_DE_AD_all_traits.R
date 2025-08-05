#Supplementary Figure 15 
rm(list=ls())
library(ggplot2)
library(data.table)
library(corrplot)
library(reshape2)
library(dplyr)

############################
####DE AD results first
################################

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
rosmap_all_results=as.data.frame(all_results)

rownames(rosmap_all_results) <- rosmap_all_results$V1
rosmap_all_results$V1 <- NULL

# 2. Pull out just the logFC columns once
rosmap_logfc_subset <- rosmap_all_results[, grep("logFC", colnames(rosmap_all_results)), drop = FALSE]

# 3. Define traits and their panel titles
traits <- c(
  "braaksc_simplified",
  "cts_mmse30",
  "ceradsc_test.txt"
)
panel_titles <- c(
  braaksc_simplified   = "Braak",
  cts_mmse30           = "MMSE30",
  "ceradsc_test.txt"   = "CERAD"
)

# 4. Standard assembly labels
assemblies <- c("NCBI36", "GRCh37", "GRCh38.p12", "GRCh38.p13", "CHM13v2.0")

# 5. Compute & store each Spearman matrix
corr_list <- list()
for (trait in traits) {
  raw_cols <- grep(trait, colnames(rosmap_logfc_subset), value = TRUE, fixed = TRUE)
  if (length(raw_cols) < 2) next
  mat <- rosmap_logfc_subset[, raw_cols, drop = FALSE]
  clean_names <- gsub("\\||\\-", "_", raw_cols)
  colnames(mat) <- clean_names
  mat <- na.omit(mat)
  
  # init with cleaned names
  res <- matrix(
    NA,
    nrow = length(clean_names),
    ncol = length(clean_names),
    dimnames = list(clean_names, clean_names)
  )
  cmb <- t(combn(clean_names, 2))
  for (i in seq_len(nrow(cmb))) {
    a <- cmb[i,1]; b <- cmb[i,2]
    sub <- mat[, c(a, b)]
    sub <- sub[!(is.na(sub[,1]) & is.na(sub[,2])), ]
    sub[is.na(sub)] <- 0
    res[a, b] <- cor(sub[,1], sub[,2], method = "spearman")
  }
  colnames(res) <- rownames(res) <- assemblies
  corr_list[[trait]] <- res
}

# 6. Plot all five in a 2×3 grid with adjusted text sizes
out_pdf <- "rosmap_cor_DE_AD_all_traits.pdf"
pdf(out_pdf, width = 14, height = 10)
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))

for (trait in traits) {
  if (! trait %in% names(corr_list)) {
    plot.new()
    title(main = panel_titles[trait], cex.main = 2.0)
    next
  }
  
  m       <- corr_list[[trait]]
  min_val <- min(m, na.rm = TRUE)
  
  corrplot(
    m,
    type        = "upper",
    method      = "color",
    col         = colorRampPalette(c("grey", "black"))(100),
    cl.lim      = c(min_val, 1),
    addgrid.col = "white",
    diag        = FALSE,
    is.corr     = FALSE,
    tl.col      = "black",
    tl.cex      = 2.0,  # axis labels
    number.cex  = 2.0,  # correlation numbers
    cl.cex      = 1.5,  # legend text
    number.digits = 2,
    cl.ratio      = 0.3
  )
  title(main = panel_titles[trait], cex.main = 2.0)
}

# 7th empty panel
plot.new()
dev.off()

