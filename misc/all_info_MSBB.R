#load libraries
library(ggplot2)
library(data.table)
library(tidyr)
library(GGally)
library(stringr)

## DE output - direct output from toptable
hg18_hg19_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg18_v30_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg19_v30_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v30_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg18_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg19_v43_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg18_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
hg19_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v30_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts
v43_T2T_DE=read.delim("DE_assembly_corrected_for_id_assembly_test.txt") #counts

### Mapping between assemblies for common gene names
mapping=readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.RDS

# Non-T2T comparison grids
pairs = c("hg18_hg19","hg18_v30","hg19_v30","v30_v43","hg18_v43","hg19_v43")
de_output = c("hg18_hg19_DE","hg18_v30_DE","hg19_v30_DE","v30_v43_DE","hg18_v43_DE","hg19_v43_DE")
combos <- c("GRCh37-NCBI36","GRCh38.p12-NCBI36","GRCh38.p12-GRCh37", "GRCh38.p13-GRCh38.p12", "GRCh38.p13-NCBI36", "GRCh38.p13-GRCh37")
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
	final_non_t2t[[x]] = de_results_subset
}

combined_df_non_T2T <- do.call(rbind, final_non_t2t)  # Combine all non-T2T

### T2T comparisons
pairs_t2t = c("hg18_T2T","hg19_T2T","v30_T2T","v43_T2T")
de_output_t2t = c("hg18_T2T_DE","hg19_T2T_DE","v30_T2T_DE","v43_T2T_DE")
combos_t2t = c("CHM13v2.0-NCBI36","CHM13v2.0-GRCh37","CHM13v2.0-GRCh38.p12","CHM13v2.0-GRCh38.p13")
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
	print("done")
	final_t2t[[x]] = de_results_subset
}

combined_df_T2T <- do.call(rbind, final_t2t)  # Combine all T2T

## Combine all information and clean for plotting
all_info_MSBB = rbind(combined_df_non_T2T, combined_df_T2T)

all_info_MSBB$direction_signif="notSignif"
all_info_MSBB$direction_signif[all_info_MSBB$sig_status=="Signif" & all_info_MSBB$direction== "New"]="New"
all_info_MSBB$direction_signif[all_info_MSBB$sig_status=="Signif" & all_info_MSBB$direction== "Old"]="Old"

all_info_MSBB$for_plot=paste0(all_info_MSBB$direction_signif,"_", all_info_MSBB$in_common )

all_info_MSBB$assembly <- factor(all_info_MSBB$assembly, levels = c("GRCh37-NCBI36", "GRCh38.p12-NCBI36", "GRCh38.p13-NCBI36", "CHM13v2.0-NCBI36", "GRCh38.p12-GRCh37", "GRCh38.p13-GRCh37", "CHM13v2.0-GRCh37", "GRCh38.p13-GRCh38.p12", "CHM13v2.0-GRCh38.p12", "CHM13v2.0-GRCh38.p13"))

saveRDS(all_info_MSBB,"all_info_MSBB.RDS" )
