

rm(list=ls())
library(data.table)

mapping=readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.R

phenos=c("braaksc_simplified","ceradsc_defvsctl","cogdx_simplified","cts_mmse30")
# assemblies=c("T2T","v43","v30","hg19","hg18")
assemblies=c("hg18","hg19","v30","v43","T2T")

assembliesComb = t(combn(assemblies,2))
assembliesComb2=assembliesComb[assembliesComb[,2]=="T2T",]
assembliesComb=assembliesComb[assembliesComb[,2]!="T2T",]
assembliesComb=rbind(assembliesComb[,2:1],assembliesComb2[,2:1])
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}),apply(assembliesComb,1,function(x){paste0(x[1],"_",x[2])}))

grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
grid$Var4[grid$Var4=="ceradsc_defvsctl"]="ceradsc_test.txt"
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies", "assembly_map")
grid=merge(grid,assembliesComb,by="assemblies")
grid$assembly_map_path=paste("path", grid$assembly1, grid$assembly2, sep = "_")
grid2=unique(grid[,c("assemblies", "assembly_map", "assembly_map_path")])

path_hg19_hg18=#path to interaction results
path_v30_hg18=#path to interaction results
path_v43_hg18=#path to interaction results
path_T2T_hg18=#path to interaction results
path_v30_hg19=#path to interaction results
path_v43_hg19=#path to interaction results
path_T2T_hg19=#path to interaction results
path_T2T_v30=#path to interaction results
path_v43_v30=#path to interaction results
path_T2T_v43=#path to interaction results


for(i in 1:nrow(grid2)){
path=get(grid2[i,"assembly_map_path"])
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
colnames(all)=paste0(colnames(all),"_",grid2[i,"assembly_map"])
colnames(all)[1]="gene_id"
value_assigned=as.character(grid2[i, "assemblies"])
assign(value_assigned,all)
print(i)
}


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
saveRDS(all_results,"results_interaction_ROSMAP_1.8.24.RDS")
