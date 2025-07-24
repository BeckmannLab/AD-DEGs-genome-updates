
#load library
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#function
format_for_gsea = function(df){
        #######input for this is a two column dataframe where the first column are gene symbols and second column DE status which is either Signif or notSignif. if you want to do gene IDs you need to switch out. output is two column dataframe where first column is gene symbols second column is signif or not signif
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

df_all = readRDS("MSBB_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/MSBB_ad_interaction_results_matrix.R

df_all2 = df_all[,c("gene_symbol", "DEAD_status", "trait")]
df_all2 <- df_all2[!grepl("CDR_defvsctl|CERJ_simplified|PlaqueMean_defvsctl|bbscore_simplified", df_all2$trait), ]
df_all2$DEAD_status = ifelse(df_all2$DEAD_status != "notSignif", "Signif", "notSignif")

colnames(df_all2) = c("ID", "DE", "trait")

#create universe and sig files
for (x in unique(df_all2$trait)){
	print(x)
	df_all_input = df_all2[df_all2$trait == x,]
	for_gsea_df = format_for_gsea(data.frame(gene_symbol = df_all_input$ID, sig = df_all_input$DE))
	for_gsea_df =unique(for_gsea_df)
	subset_df2 = for_gsea_df$gene_ID
	write.table(subset_df2,paste0(x,"_msigdb_singlecell_universe.txt"), row.names = FALSE)
	sig_subset_df = for_gsea_df[which(for_gsea_df$DE == "Signif"),]
	sig_subset_df2 =sig_subset_df$gene_ID
	write.table(sig_subset_df2,paste0(x,"_msigdb_singlecell_sig.txt"),row.names = FALSE)
}


##c2 
for (x in unique(df_all2$trait)){
	print(x)
	var1=paste0(x,"_msigdb_singlecell_sig.txt")
	var2=paste0(x,"_msigdb_singlecell_universe.txt")
	var3="c2.all.v2023.2.Hs.symbols.gmt"  #from msigdb website C2 curated gene sets  https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
	var4="/c2/"
	var5=paste0(x, "_single_cell_msigdb.txt")
	system(paste("python enrich_function.py",var1, var2, var3, var4, var5)) ##see AD-DEGs-genome-updates/misc/enrich_function.py
}

folder <- "/c2/"

# get full paths to all .csv files there
files <- list.files(path = folder, pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(files, function(f) {
  # read in (adjust read.* as needed)
  df <- fread(f)
  
  # grab basename, then capture first two pieces before '_'
  trait <- sub("^([^_]+_[^_]+).*", "\\1", basename(f))
  
  # add as a new column
  df$trait <- trait
  
  df
})

# combine into one data.frame
df <- do.call(rbind, data_list)
df2_c2 = as.data.frame(df)
df2_c2$trait[df2_c2$trait == "PlaqueMean_single"] <- "PlaqueMean"

saveRDS(df2_c2, "c2_msigdb_msbb.RDS")

##hallmark 
for (x in unique(df_all2$trait)){
	print(x)
	var1=paste0(x,"_msigdb_singlecell_sig.txt")
	var2=paste0(x,"_msigdb_singlecell_universe.txt")
	var3="h.all.v2025.1.Hs.symbols.gmt" #from msigdb website hallmark gene sets https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
	var4="/hallmark"
	var5=paste0(x, "_single_cell_msigdb.txt")
	system(paste("python enrich_function.py",var1, var2, var3, var4, var5)) ##see AD-DEGs-genome-updates/misc/enrich_function.py
}

folder <- "/hallmark"

# get full paths to all .csv files there
files <- list.files(path = folder, pattern = "\\.txt$", full.names = TRUE)

data_list <- lapply(files, function(f) {
  # read in (adjust read.* as needed)
  df <- fread(f)
  
  # grab basename, then capture first two pieces before '_'
  trait <- sub("^([^_]+_[^_]+).*", "\\1", basename(f))
  
  # add as a new column
  df$trait <- trait
  
  df
})

# combine into one data.frame
df <- do.call(rbind, data_list)
df2_hallmark = as.data.frame(df)
df2_hallmark$trait[df2_hallmark$trait == "PlaqueMean_single"] <- "PlaqueMean"
all = rbind(df2_hallmark, df2_c2)

#adjust p value
all$correct_adj_p_value = p.adjust(all$'P-value', method = "BH")
subset_df = all[all$correct_adj_p_value <= 0.05,]
subset_df = subset_df[order(subset_df$correct_adj_p_value, decreasing = FALSE),]

saveRDS(all, "c2_hallmark_msigdb_msbb.RDS") #Supplementary Table 2: MSBB msigdb across ref pairs
