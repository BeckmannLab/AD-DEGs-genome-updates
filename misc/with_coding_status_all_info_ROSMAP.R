## Clear workspace and load required libraries
rm(list=ls())

#load libraries
library(ggplot2)
library(data.table)
library(tidyr)
library(GGally)
library(stringr)

## DE output - direct output from toptable
hg18_hg19_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_v30_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_v30_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v30_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v30_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v43_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt")

# GTF files
hg18_gtf=fread("Homo_sapiens.NCBI36.53.gtf")        # Load hg18 GTF annotation
hg19_gtf=fread("Homo_sapiens.GRCh37.70.gtf")        # Load hg19 GTF annotation
v30_gtf=fread("gencode.v30.primary_assembly.annotation.gtf")     # Load GRCh38 v30 GTF
v43_gtf=fread("gencode.v43.primary_assembly.annotation.gtf", data.table= FALSE)  # Load GRCh38 v43 GTF
t2t_gtf=fread("Homo_sapiens_GCA_009914755.4_2022_07_genes.gtf",  data.table = FALSE)  # Load T2T GTF

### Mapping between assemblies for common gene names
mapping=readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.RDS

## Pre-formatting

# hg18
colnames(hg18_gtf)=c("genomic region","type","region","start","end","dot","strand","not_sure","V9")
split_columns <- data.frame(do.call("rbind", strsplit(as.character(hg18_gtf$V9), ";",fixed = TRUE))) 
result_df <- cbind(hg18_gtf, split_columns)
sub_hg18=result_df[,c("X1", "type")]   
hg18_df=data.frame(do.call("rbind", strsplit(as.character(sub_hg18$X1), " ",fixed = TRUE)))  
result_df <- cbind(hg18_df,sub_hg18$type)
result_df$X2 <- gsub('"', '', result_df$X2)   
hg18_final=result_df[,2:3]
colnames(hg18_final)=c("gene_id", "status")  
hg18_final$status[hg18_final$status != "protein_coding"] <- "non_coding"  
hg18_final=unique(hg18_final)  

# hg19
colnames(hg19_gtf)=c("genomic region","type","region","start","end","dot","strand","not_sure","V9")
split_columns <- data.frame(do.call("rbind", strsplit(as.character(hg19_gtf$V9), ";",fixed = TRUE)))
result_df <- cbind(hg19_gtf, split_columns)
result_df=result_df[which(result_df$region=="exon"),]  
sub_hg19=result_df[,c("X1", "X5")] 
sub_hg19$status <- str_replace_all(str_extract_all(sub_hg19$X5, '"([^"]+)"'), '"', '') 
sub_hg19$gene_id <- str_replace_all(str_extract_all(sub_hg19$X1, '"([^"]+)"'), '"', '')  
dag19_df=sub_hg19[,c("gene_id", "status")] 
hg19_df=unique(hg19_df)
hg19_df$status[hg19_df$status != "protein_coding"] <- "non_coding"  
hg19_final=unique(hg19_df)  

# v30
colnames(v30_gtf)=c("genomic region","type","region","start","end","dot","strand","not_sure","V9")
result_df=v30_gtf[which(v30_gtf$region=="exon"),] 
result <- sub(".*gene_type\\s(.*?);.*", "\\1", result_df$V9)  
result2 <- sub(".*gene_id\\s(.*?);.*", "\\1", result_df$V9)  
result_df2 <- cbind(result2, result)
result_df2 <- as.data.frame(result_df2)
result_df2$result <- gsub('"', '', result_df2$result) 
result_df2$result2 <- gsub('"', '', result_df2$result2)  
colnames(result_df2)=c("gene_id", "status")
result_df2$status[result_df2$status != "protein_coding"] <- "non_coding"
v30_final=result_df2
v30_final=unique(v30_final)
v30_final$Geneid_noVersion=unlist(lapply(strsplit(as.character(v30_final$gene_id),".",fixed=TRUE),function(x){ 
		if(sum(grepl("_", x))==FALSE){
			x[1]
		}else{
			paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
			}
		}))
v30_final=v30_final[,c(3,2)]
colnames(v30_final)=c("gene_id", "status") 

# v43
colnames(v43_gtf)=c("genomic region","type","region","start","end","dot","strand","not_sure","V9")
v43=v43_gtf[which(v43_gtf$region=="exon"),]
result <- sub(".*gene_type\\s(.*?);.*", "\\1", v43$V9)
result2 <- sub(".*gene_id\\s(.*?);.*", "\\1", v43$V9)
result_df <- cbind(result2, result)
result_df <- as.data.frame(result_df)
result_df$result <- gsub('"', '', result_df$result)
result_df$result2 <- gsub('"', '', result_df$result2)
colnames(result_df)=c("gene_id", "status")
result_df$status[result_df$status != "protein_coding"] <- "non_coding"
v43_final=result_df
v43_final=unique(v43_final)
v43_final$Geneid_noVersion=unlist(lapply(strsplit(as.character(v43_final$gene_id),".",fixed=TRUE),function(x){
		if(sum(grepl("_", x))==FALSE){
			x[1]
		}else{
			paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
			}
		}))
v43_final=v43_final[,c(3,2)]
colnames(v43_final)=c("gene_id", "status")

# t2t
colnames(t2t_gtf)=c("genomic region","type","region","start","end","dot","strand","not_sure","V9")
t2t=t2t_gtf[which(t2t_gtf$region=="exon"),]
result <- sub(".*gene_biotype\\s(.*?);.*", "\\1", t2t$V9)
result2 <- sub(".*gene_id\\s(.*?);.*", "\\1", t2t$V9)
result_df <- cbind(result2, result)
result_df <- as.data.frame(result_df)
result_df$result <- gsub('"', '', result_df$result)
result_df$result2 <- gsub('"', '', result_df$result2)
colnames(result_df)=c("gene_id", "status")
result_df$status[result_df$status != "protein_coding"] <- "non_coding"
t2t_final=result_df
t2t_final=unique(t2t_final)  

# Non-T2T comparison grids
pairs = c("hg18_hg19","hg18_v30","hg19_v30","v30_v43","hg18_v43","hg19_v43")
de_output = c("hg18_hg19_DE","hg18_v30_DE","hg19_v30_DE","v30_v43_DE","hg18_v43_DE","hg19_v43_DE")
combos <- c("GRCh37-NCBI36","GRCh38.12-NCBI36","GRCh38.12-GRCh37", "GRCh38.13-GRCh38.12", "GRCh38.13-NCBI36", "GRCh38.13-GRCh37")
gtf_files1 = c("hg19_final", "v30_final", "v30_final", "v43_final", "v43_final", "v43_final")
gtf_files2 = c("hg18_final", "hg18_final","hg19_final", "v30_final", "hg18_final", "hg19_final")
grid = as.data.frame(cbind(de_output,pairs,combos,gtf_files1,gtf_files2))

# Loop for non-T2T
final_non_t2t <- vector("list", nrow(grid))
for (x in seq_len(nrow(grid))){
	print(x)  
	# DE
	de_results = get(grid[x,1])

	# in common mapping
	reference_common=mapping[grep(grid[x,2], mapping$combo),]
	reference_common <- na.omit(reference_common)
  	de_reference_common=de_results[de_results$X %in% reference_common$common_name,]
	print("common")

	# not in common mapping
	reference_not_common=mapping[grep(grid[x,2], mapping$combo),]
	reference_not_common <- reference_not_common[!complete.cases(reference_not_common), ]
  	de_reference_not_common=de_results[de_results$X %in% reference_not_common$common_name,]
	print("not in common")

  	# Annotate DE results
  	de_results$sig_status="notSignif"
	de_results$sig_status[de_results$adj.P.Val < 0.05 ]="Signif"
	de_results$direction[de_results$logFC < 0 ]="Old"
	de_results$direction[de_results$logFC > 0 ]="New"
	de_results$in_common=TRUE
	print("DE clean")
  	de_results$in_common[de_results$X %in% reference_not_common$common_name]=FALSE
  	de_results$assembly=grid[x,3]  
	colnames(de_results)[1]="gene_id"
	de_results_subset = de_results[,c("gene_id", "sig_status", "direction", "in_common", "assembly")]
	print("combine")

	# Fetch coding status tables
	ref_final1 = get(grid[x,4])
    print("ref1")

	ref_final2 = get(grid[x,5])
    print("ref2")

    # Merge coding status and DE
    both_ref_status=merge(ref_final1, ref_final2, all = T)
	ref_new=merge(both_ref_status, de_results_subset, by = "gene_id", all = T)
	print("done")
	final_non_t2t[[x]] = ref_new
}

combined_df_non_T2T <- do.call(rbind, final_non_t2t)  # Combine all non-T2T

### T2T comparisons
pairs_t2t = c("hg18_T2T","hg19_T2T","v30_T2T","v43_T2T")
de_output_t2t = c("hg18_T2T_DE","hg19_T2T_DE","v30_T2T_DE","v43_T2T_DE")
combos_t2t = c("CHM13v2.0-NCBI36","CHM13v2.0-GRCh37","CHM13v2.0-GRCh38.12","CHM13v2.0-GRCh38.13")
gtf_files1_t2t = c("t2t_final", "t2t_final", "t2t_final", "t2t_final")
gtf_files2_t2t = c("hg18_final", "hg19_final", "v30_final", "v43_final")
grid_t2t = as.data.frame(cbind(de_output_t2t,pairs_t2t,combos_t2t,gtf_files1_t2t,gtf_files2_t2t))

# Loop for T2T
final_t2t <- vector("list", nrow(grid_t2t))
for (x in seq_len(nrow(grid_t2t))){
	print(x)
	# DE
	de_results = get(grid_t2t[x,1])

	# common vs not common mapping
	reference_common=mapping[grep(grid_t2t[x,2], mapping$combo),]
	reference_common <- na.omit(reference_common)
  	de_reference_common=de_results[de_results$X %in% reference_common$common_name,]
	print("common")

	reference_not_common=mapping[grep(grid_t2t[x,2], mapping$combo),]
	reference_not_common <- reference_not_common[!complete.cases(reference_not_common), ]
  	de_reference_not_common=de_results[de_results$X %in% reference_not_common$common_name,]
	print("not in common")

  	# Annotate DE
  	de_results$sig_status="notSignif"
	de_results$sig_status[de_results$adj.P.Val < 0.05 ]="Signif"
	de_results$direction[de_results$logFC < 0 ]="Old"
	de_results$direction[de_results$logFC > 0 ]="New"
	de_results$in_common=TRUE
	print("DE clean")
  	de_results$in_common[de_results$X %in% reference_not_common$common_name]=FALSE
  	de_results$assembly=grid_t2t[x,3]
	colnames(de_results)[1]="gene_id"
	de_results_subset = de_results[,c("gene_id", "sig_status", "direction", "in_common", "assembly")]
	print("combine")

	# Merge T2T gene status via mapping
	t2t_map=mapping[which(mapping$combo==grid_t2t[x,2]),]
	t2t_final2=merge(t2t_final, t2t_map, by.x="gene_id", by.y="assembly_2_orginal_gene_name", all = T)
	t2t_both=t2t_final2[,c("common_name", "status")]
	t2t_both=unique(t2t_both)
	colnames(t2t_both)=c("gene_id", "status")
	both_status=merge(get(grid_t2t[x,5]), t2t_both, all = T)  # Merge with second GTF
	both_status_new=merge(both_status, de_results_subset, by="gene_id", all = T)

	print("done")
	final_t2t[[x]] = both_status_new
}

combined_df_T2T <- do.call(rbind, final_t2t)  # Combine all T2T

## Combine all information and clean for plotting
all_info_ROSMAP = rbind(combined_df_non_T2T, combined_df_T2T)
df_clean <- all_info_ROSMAP[complete.cases(all_info_ROSMAP$in_common, all_info_ROSMAP$assembly, all_info_ROSMAP$direction, all_info_ROSMAP$sig_status), ]
df_clean2 <- df_clean[complete.cases(df_clean$status), ]
all_info_ROSMAP=df_clean2
# Create combined direction_signif field
all_info_ROSMAP$direction_signif="notSignif"
all_info_ROSMAP$direction_signif[all_info_ROSMAP$sig_status=="Signif" & all_info_ROSMAP$direction== "New"]="New"
all_info_ROSMAP$direction_signif[all_info_ROSMAP$sig_status=="Signif" & all_info_ROSMAP$direction== "Old"]="Old"
# Prepare factor levels and plotting
all_info_ROSMAP$assembly <- factor(all_info_ROSMAP$assembly, levels = c("GRCh37-NCBI36", "GRCh38.12-NCBI36", "GRCh38.13-NCBI36", "CHM13v2.0-NCBI36", "GRCh38.12-GRCh37", "GRCh38.13-GRCh37", "CHM13v2.0-GRCh37", "GRCh38.13-GRCh38.12", "CHM13v2.0-GRCh38.12", "CHM13v2.0-GRCh38.13"))

## Save final data
saveRDS(all_info_ROSMAP,"with_coding_status_all_info_ROSMAP.RDS" )
