library(dplyr)
library(rex)
library(tidyr)
library(stringr)

df <- readRDS("MSBB_interaction_andAD_DE_with_symbols_4.10.24.RDS") # see AD_DEGs-genome-updates/misc/MSBB_interaction_andAD_DE_with_symbols_4.10.24.R

assemblies=c("hg18","hg19","v30","v43","T2T")
mart=read.delim("mart_export.txt") #see AD-DEGs-genome-updates/files/mart_export.txt

colnames(mart)[1]="gene_id"
colnames(mart)[4]="gene_symbol"
colnames(df) = gsub("PlaqueMean_test.txt","PlaqueMean", colnames(df))

########################
#######################

# Specify phenotypes or tests to process
phenos=c("CDR_simplified", "CERJ_defvsctl","PlaqueMean_test.txt")

# List of genome assemblies
assemblies=c("hg18","hg19","v30","v43","T2T")

## Generate all unique pairs of assemblies
assembliesComb = t(combn(assemblies,2))
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}))

# Create a grid of assembly pairs and phenotypes
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")
grid$combo= paste0(grid$assembly2, "_", grid$assembly1)

##pull logFC
df_subset_logFC_DEAD = df %>% select(gene_id, contains("DEAD_logFC"))
df_subset_logFC_interaction = df %>% select(gene_id, contains("interaction_logFC"))

##pull status 
df_subset_DEAD_status = df %>% select(gene_id, matches(rex(start, "DEAD_", one_or_more(none_of("_")), "DE_", one_or_more(any), end)))
df_subset_interaction_status = df %>% select(gene_id, matches(rex(start, "interaction_", one_or_more(none_of("_")), "DE_", one_or_more(any), end)))

##pull p values 
df_subset_DEAD_pval = df %>% 
  select(gene_id, contains("DEAD_adj.P.Val"))

df_subset_interaction_pval = df %>% 
  select(gene_id, contains("interaction_adj.P.Val"))

#format
colnames(df_subset_logFC_DEAD) = gsub("\\.", "_",colnames(df_subset_logFC_DEAD))
colnames(df_subset_DEAD_pval) = gsub("\\.", "_",colnames(df_subset_DEAD_pval))
colnames(df_subset_interaction_pval) = gsub("\\.", "_",colnames(df_subset_interaction_pval))
colnames(df_subset_logFC_interaction) = gsub("\\.", "_",colnames(df_subset_logFC_interaction))

##################
#DEAD edits
##################

##edit DEAD ones
df_subset_logFC_DEAD_long = pivot_longer(
	df_subset_logFC_DEAD, 
	contains("DEAD_logFC"),
	names_to="test",
	values_to="DEAD_logFC",
	names_prefix=rex("DEAD_logFC|")
) %>% 
	mutate(
		test = str_replace(test, rex("PATH.Dx_test.txt"), "PATH_Dx_test")
	) %>%
	mutate(
		test = str_replace(test, rex(capture(anything), "_"), "\\1.")
	) %>%
	separate(
		test,
		into=c("trait", "assembly"),
		sep =rex(".")
	)

df_subset_pval_DEAD_long = pivot_longer(
	df_subset_DEAD_pval, 
	contains("DEAD_adj_P_Val"),
	names_to="test",
	values_to="DEAD_adj_P_Val",
	names_prefix=rex("DEAD_adj_P_Val|")
) %>% 
	mutate(
		test = str_replace(test, rex("PATH.Dx_test.txt"), "PATH_Dx_test")

	) %>%
	mutate(
		test = str_replace(test, rex(capture(anything), "_"), "\\1.")
	) %>%
	separate(
		test,
		into=c("trait", "assembly"),
		sep =rex(".")
	)

df_subset_logFC_DEAD_long2 = merge(df_subset_logFC_DEAD_long,df_subset_pval_DEAD_long, by = c("gene_id", "trait", "assembly") )


df_subset_logFC_DEAD_diff=inner_join(
	df_subset_logFC_DEAD_long2, 
	df_subset_logFC_DEAD_long2, 
	by = c("gene_id", "trait"), 
	relationship = "many-to-many", 
	suffix = c(".older", ".newer")
) %>%
	mutate(
		assembly_comparison=paste0(assembly.newer, "_", assembly.older),
		DE_AD_difference=abs(DEAD_logFC.newer-DEAD_logFC.older),
		DE_AD_sign_same=sign(DEAD_logFC.newer) * sign(DEAD_logFC.older) >= 0,
	) %>%
	filter(
		assembly_comparison %in% grid$combo
		)

df_subset_DEAD_status_long= pivot_longer(
	df_subset_DEAD_status, 
	starts_with("DEAD_"), 
	names_to = "test", 
	values_to="DEAD_status", 
	names_prefix="DEAD_"
)%>%
	separate(
		test, 
		into=c("assembly_comparison", "trait"),
		sep =rex("DE_")
	) %>%
	mutate(
		assembly_comparison=str_replace(
			assembly_comparison,
			rex(start, capture(or(assemblies)), capture(or(assemblies)), end), 
			"\\1_\\2"
		)
	)

df_subset_logFC_DEAD_long2 = merge(df_subset_logFC_DEAD_diff,df_subset_pval_DEAD_long, by = c("gene_id", "trait") )

##edit interaction ones
df_subset_logFC_interaction_long = pivot_longer(
	df_subset_logFC_interaction, 
	contains("interaction_logFC"),
	names_to="test",
	values_to="interaction_logFC",
	names_prefix=rex("interaction_logFC|")
)%>%
	mutate(
    test = str_replace(test, rex(capture(anything), "_", capture(anything), "_"), "\\1.\\2_")
  ) %>%
	separate(
		test,
		into=c("trait", "assembly_comparison"),
		sep =rex(".")
	)

df_subset_interaction_status_long= pivot_longer(
	df_subset_interaction_status, 
	starts_with("interaction_"), 
	names_to = "test", 
	values_to="interaction_status", 
	names_prefix="interaction_"
)%>%
	separate(
		test, 
		into=c("assembly_comparison", "trait"),
		sep =rex("DE_")
	) %>%
	mutate(
		assembly_comparison=str_replace(
			assembly_comparison,
			rex(start, capture(or(assemblies)), capture(or(assemblies)), end), 
			"\\1_\\2"
		),
		interaction_status=if_else(interaction_status=="notSignif", "notSignif", "Signif")
	)

df_subset_pval_interaction_long = pivot_longer(
	df_subset_interaction_pval, 
	contains("interaction_adj_P_Val"),
	names_to="test",
	values_to="interaction_adj_P_Val",
	names_prefix=rex("interaction_adj_P_Val|")
)%>%
	mutate(
    test = str_replace(test, rex(capture(anything), "_", capture(anything), "_"), "\\1.\\2_")
  ) %>%
	separate(
		test,
		into=c("trait", "assembly_comparison"),
		sep =rex(".")
	)

##combine them

all = merge(df_subset_pval_interaction_long, df_subset_pval_DEAD_long, by = c("gene_id", "trait"))


df_subset_interaction_status_long2 = merge(df_subset_interaction_status_long,df_subset_pval_interaction_long, by = c("gene_id", "trait", "assembly_comparison") )


df_all=df_subset_logFC_DEAD_diff %>% 
	inner_join(
		df_subset_DEAD_status_long, 
		relationship = "one-to-one",
		by = c("gene_id", "trait", "assembly_comparison")
	) %>%
	inner_join(
		df_subset_interaction_status_long2, 
		relationship = "one-to-one",
		by = c("gene_id", "trait", "assembly_comparison")
	) %>%
	inner_join(
		mart, 
		by = "gene_id"
	)%>%
	arrange(assembly_comparison, trait, gene_id)

df_all=as.data.frame(df_all)

saveRDS(df_all, "MSBB_ad_interaction_results_matrix.RDS")
