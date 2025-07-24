## Clear the global environment
rm(list = ls())

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
phenos=c("CDR_simplified", "CERJ_defvsctl","PlaqueMean_test.txt")

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
  # Read and merge each DE file
  for(DE_file in DE_files){
      tmp=fread(paste0(path,DE_file))
      phenotype=unlist(lapply(strsplit(DE_file,"_"),function(x){paste0(x[11],"_",x[12])}))
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
saveRDS(all_results, "MSBB_with_signif_categories_results_AD_all_resultsDE_4.10.24.RDS")

#######################
###interaction ones
########################
library(data.table)

mapping=readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.R

# List of genome assemblies
assemblies=c("hg18","hg19","v30","v43","T2T")

# Specify phenotypes or tests to process
phenos <- c("CDR_simplified", "CERJ_defvsctl", "PlaqueMean_test.txt")

# Generate all unique pairs of assemblies
assembliesComb = t(combn(assemblies,2))
assembliesComb2=assembliesComb[assembliesComb[,2]=="T2T",]
assembliesComb=assembliesComb[assembliesComb[,2]!="T2T",]
assembliesComb=rbind(assembliesComb[,2:1],assembliesComb2[,2:1])
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}),apply(assembliesComb,1,function(x){paste0(x[1],"_",x[2])}))

# Create a grid of assembly pairs and phenotypes
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
grid$Var4[grid$Var4=="ceradsc_defvsctl"]="ceradsc_test.txt"
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies", "assembly_map")
grid=merge(grid,assembliesComb,by="assemblies")
grid$assembly_map_path=paste("path", grid$assembly1, grid$assembly2, sep = "_")
grid2=unique(grid[,c("assemblies", "assembly_map", "assembly_map_path")])

# Define file paths for differential expression results for each genome assembly
path_hg19_hg18= #path to interaction results
path_v30_hg18= #path to interaction results
path_v43_hg18= #path to interaction results
path_T2T_hg18= #path to interaction results
path_v30_hg19= #path to interaction results
path_v43_hg19= #path to interaction results
path_T2T_hg19= #path to interaction results
path_T2T_v30= #path to interaction results
path_v43_v30= #path to interaction results
path_T2T_v43= #path to interaction results

# Loop through each assembly to read and merge DE files
for(i in 1:nrow(grid2)){
path=get(grid2[i,"assembly_map_path"])
DE_files=list.files(path,pattern="^DE")
all=data.table()
    for(DE_file in DE_files){
        tmp=fread(paste0(path,DE_file))
        phenotype=unlist(lapply(strsplit(DE_file,"_"),function(x){paste0(x[11],"_",x[12])}))
        colnames(tmp)[2:ncol(tmp)]=paste0(colnames(tmp)[2:ncol(tmp)],"|",phenotype)
        if(nrow(all)<1){
            all=tmp
        }else{
            all=merge(all,tmp,by="V1",all=TRUE)
        }}
colnames(all)=paste0(colnames(all),"_",grid2[i,"assembly_map"])
colnames(all)[1]="gene_id"
value_assigned=as.character(grid2[i, "assemblies"])
assign(value_assigned,all)
}

# Merge results from all assemblies into one data frame
all_results=merge(hg19hg18,v30hg18,by="gene_id",all=TRUE)
all_results=merge(all_results,v43hg18,by="gene_id",all=TRUE)
all_results=merge(all_results,v30hg19,by="gene_id",all=TRUE)
all_results=merge(all_results,v43hg19,by="gene_id",all=TRUE)
all_results=merge(all_results,v43v30,by="gene_id",all=TRUE)
all_results=merge(all_results,T2Thg18,by="gene_id",all=TRUE)
all_results=merge(all_results,T2Thg19,by="gene_id",all=TRUE)
all_results=merge(all_results,T2Tv30,by="gene_id",all=TRUE)
all_results=merge(all_results,T2Tv43,by="gene_id",all=TRUE)
all_results=as.data.frame(all_results)

# Reinitialize combinations for significance categorization
assembliesComb = t(combn(assemblies,2))
assembliesComb2=assembliesComb[assembliesComb[,2]=="T2T",]
assembliesComb=assembliesComb[assembliesComb[,2]!="T2T",]
assembliesComb=rbind(assembliesComb[,2:1],assembliesComb2[,2:1])
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}),apply(assembliesComb,1,function(x){paste0(x[2],"_",x[1])}))

# Recreate grid for significance loop
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies", "assembly_map")
grid=merge(grid,assembliesComb,by="assemblies")

# Copy merged results for analysis
interaction_msbb = all_results

# Loop to categorize significance across assembly pairs for each phenotype
for(i in 1:nrow(grid)){
    outCol=paste0(grid[i,"assembly1"],grid[i,"assembly2"],"DE_",grid[i,"pheno"])
    colOfInterest1=paste0("adj.P.Val|",grid[i,"pheno"],"_",grid[i,"assembly1"],"_",grid[i,"assembly2"])
    colOfInterest2=paste0("logFC|",grid[i,"pheno"],"_",grid[i,"assembly1"],"_",grid[i,"assembly2"])
    interaction_msbb[,outCol]=""
    interaction_msbb[!is.na(interaction_msbb[,colOfInterest1]),outCol]="notSignif"
    # interaction_msbb[which(interaction_msbb[,colOfInterest1]<0.05 & interaction_msbb[,colOfInterest2]>0.05),outCol]="bothSignif"
    interaction_msbb[which(interaction_msbb[,colOfInterest1]<=0.05 & interaction_msbb[,colOfInterest2]> 0),outCol]=paste0(grid[i,"assembly1"],"","Signif")
    interaction_msbb[which(interaction_msbb[,colOfInterest1]<=0.05 & interaction_msbb[,colOfInterest2]<= 0),outCol]=paste0(grid[i,"assembly2"],"Signif")
    interaction_msbb[interaction_msbb[,outCol]=="",outCol]=NA

    interaction_msbb[,outCol]=factor(interaction_msbb[,outCol],levels=c(paste0(grid[i,"assembly1"],"Signif"),paste0(grid[i,"assembly2"],"Signif"),"notSignif"))

    ##subset for ones in common
    mapping_subset=mapping[mapping$combo==grid$assembly_map[i],] ##pick the right assembly combo 

    mapping_subset2 = mapping_subset[rowSums(is.na(mapping_subset)) > 0, ] ##which ones have NA = not in common
    interaction_msbb[interaction_msbb$gene_id %in% mapping_subset2$common_name,outCol]=NA ##NA means not in common

}

# save results
saveRDS(interaction_msbb, "with_signif_categories_results_interaction_MSBB_4.7.24.RDS")

##Load them in to combine interaction and DE AD results
AD_DE=readRDS("MSBB_with_signif_categories_results_AD_all_resultsDE_4.10.24.RDS")
interaction_DE=readRDS("with_signif_categories_results_interaction_MSBB_4.7.24.RDS")

###merge them
colnames(AD_DE)=paste("DEAD", colnames(AD_DE), sep = "_")
colnames(interaction_DE)=paste("interaction", colnames(interaction_DE), sep = "_")
colnames(AD_DE)[1]="gene_id"
colnames(interaction_DE)[1]="gene_id"
interaction_and_AD=merge(AD_DE, interaction_DE, by = "gene_id", all = T)

##add symbol
mart=read.delim("mart_export.txt") #see AD-DEGs-genome-updates/files/mart_export.txt
colnames(mart)[1]="gene_id"
colnames(mart)[4]="gene_symbol"

interaction_and_AD_with_symbols=merge(interaction_and_AD, mart, by = "gene_id")

saveRDS(interaction_and_AD_with_symbols, "MSBB_interaction_andAD_DE_with_symbols_4.10.24.RDS")

