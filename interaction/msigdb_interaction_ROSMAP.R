#load library
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


format_for_gsea = function(df){
	#######input_rosmap for this is a two column dataframe where the first column are gene symbols and second column DE status which is either Signif or notSignif. if you want to do gene IDs you need to switch out. output is two column dataframe where first column is gene symbols second column is signif or not signif
	library(assertthat)
	new_df = data.frame(gene_symbol= unique(df$gene_symbol), sig = "notSignif") ##select the columns of interest and only keep unique genes
	unique_DEGs = unique(df$gene_symbol[which(df$sig == "Signif")])
	assert_that(identical(new_df$gene_symbol[match(unique_DEGs,new_df$gene_symbol)],unique_DEGs)) #make sure they match
	
	new_df$sig[match(unique_DEGs,new_df$gene_symbol)] = "Signif" #ones that are unique become unique in the df

	assert_that(identical(sort(new_df$gene_symbol[new_df$sig=="Signif"]),sort(unique_DEGs)))
	colnames(new_df) = c("gene_ID","DE")

	##remove empty rows 
	new_df[which(new_df$gene_ID == ""),] = NA #remove symbols that dont match
	new_df =na.omit(new_df)
	new_df
}


###Load in 

df_all_rosmap = readRDS("rosmap_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/rosmap_ad_interaction_results_matrix.R
common = df_all_rosmap[which(!is.na(df_all_rosmap$DE_AD_sign_same)),]
common2= common[common$DEAD_status != "notSignif",]



##make the files that are input_rosmap
for (i in unique(common$assembly_comparison)){
	print(i)
	df_sig3 = common2[which(common2$assembly_comparison == i),]
	for_gsea_df = format_for_gsea(data.frame(gene_symbol = df_sig3$gene_symbol, sig = df_sig3$interaction_status))
	for_gsea_df =unique(for_gsea_df)
	subset_df2 = for_gsea_df$gene_ID
	subset_df2 = as.data.frame(subset_df2)
	write.table(subset_df2, paste0("rosmap_common_",i, "_background_5.13.24.txt"), row.names = F, sep = ",")
	sig_subset_df = for_gsea_df[which(for_gsea_df$DE == "Signif"),]
	sig_subset_df2 =sig_subset_df$gene_ID
	sig_subset_df2 = as.data.frame(sig_subset_df2)
	write.table(sig_subset_df2, paste0("rosmap_common_only_sig_",i, "_5.13.24.txt"), row.names = F, sep = ",")
}


##run enrichment for c1
for (x in unique(common$trait)){
	for (i in unique(common$assembly_comparison)){
	var1=paste0("rosmap_common_only_sig_",i, "_5.13.24.txt")
	var2=paste0("rosmap_common_",i, "_background_5.13.24.txt")
	var3="c1.all.v2023.2.Hs.symbols.gmt" #from msigdb website C1 curated gene sets  https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
	var4="/c1/"
	var5=paste0(x,"_", i, "msigdb")
	system(paste("python enrich_function.py",var1, var2, var3, var4, var5)) ##see AD-DEGs-genome-updates/misc/enrich_function.py
}
}

##pull

library(data.table)

path="/c1/"
output_files=list.files("/c1/")
all=data.table()
    for(x in output_files){
        tmp=fread(paste0(path,x))
        tmp$trait = x
        all=rbind(all,tmp)
        }

all = all[order(all$`Adjusted P-value`, decreasing = F),]

saveRDS(all,"rosmap_c1_msigdb.RDS")
