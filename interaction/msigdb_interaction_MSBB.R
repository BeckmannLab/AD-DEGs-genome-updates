#load library
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


format_for_gsea = function(df){
	#######input_msbb for this is a two column dataframe where the first column are gene symbols and second column DE status which is either Signif or notSignif. if you want to do gene IDs you need to switch out. output is two column dataframe where first column is gene symbols second column is signif or not signif
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

df_all_msbb = readRDS("MSBB_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/MSBB_ad_interaction_results_matrix.R
common = df_all_msbb[which(!is.na(df_all_msbb$DE_AD_sign_same)),]

####setup

traits = c("CDR_simplified", "CERJ_decsctl", "PlaqueMean")
pairs = c("hg19_hg18","v30_hg18","v30_hg18", "v43_hg18","T2T_hg18","v30_hg19", "v43_hg19", "T2T_hg19", "v43_v30", "T2T_v30", "T2T_v43")

all = expand.grid(traits, pairs)
all$filename = paste0(traits, "_",pairs, "msigDB_chromsome", ".txt")

##make the files that are input_msbb
for (i in 1:nrow(all)){
	df_sig2 = common[which(common$trait == all[i,1]),]
	df_sig3 = df_sig2[which(df_sig2$assembly_comparison == all[i,2]),]
	for_gsea_df = format_for_gsea(data.frame(gene_symbol = df_sig3$gene_symbol, sig = df_sig3$interaction_status))
	back = as.data.frame(for_gsea_df$gene_ID)
	write.table(back, paste0("/input_msbb/msbb_common_",all[i,1],"_", all[i,2], "_background_5.13.24.txt"), row.names = T, sep = ",")
	df = as.data.frame(for_gsea_df$gene_ID[for_gsea_df$DE == "Signif"])
	write.table(df, paste0("/input_msbb/","msbb_common_only_sig","_",all[i,1],"_",all[i,2], "_5.13.24.txt"), row.names = T, sep = ",")
}

##set up
gene_list_files=list.files("/input_msbb",pattern="msbb_common_only")
universe_files=list.files("/input_msbb",pattern="background")
all2 = cbind(gene_list_files, universe_files)
all2 = as.data.frame(all2)
all2$file_name = gsub("msbb_common_only_sig_", "",all2$gene_list_files)
all2$file_name = gsub("_5.13.24.txt", "",all2$file_name)

##run enrichment
for (i in 1:nrow(all2)){
	var1=paste0("/input_msbb/",all2[i,1])
	var2=paste0("/input_msbb/",all2[i,2])
	var3="c1.all.v2023.2.Hs.symbols.gmt" #from msigdb website C1 curated gene sets  https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
	var4="/output_msbb"
	var5=all2[i,3]
	system(paste("python enrich_function.py",var1, var2, var3, var4, var5)) ##see AD-DEGs-genome-updates/misc/enrich_function.py
}

library(data.table)

path="/output_msbb/"
output_files=list.files("/output_msbb/")
all=data.table()
    for(x in output_files){
        tmp=fread(paste0(path,x))
        tmp$trait = x
        all=rbind(all,tmp)
        }

all = all[order(all$`Adjusted P-value`, decreasing = F),]

saveRDS(all,"msbb_final_chromsome.RDS")
