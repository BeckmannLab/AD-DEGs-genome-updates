rm(list=ls()) #this lists all of the objects in your workspace and removes objects in the workspace. clears all the objects from teh workspace to start with a clean environment
library(data.table)
library(limma)
library(edgeR)
library(variancePartition)
library(ggplot2)
library(BiocParallel)
library(assertthat)
library(tidyverse)
library(dplyr)

###functions
	#function BatchtoolsParam allow job customization 
	MyBatchtoolsLSFParam <- function(..., registryargs = batchtoolsRegistryargs(),
	                                  conf.file = "/.batchtools.conf.R") {
	    if (missing(conf.file)) {
	        return(BatchtoolsParam(..., registryargs = registryargs))
	    }
	    registryargs$conf.file <- conf.file
	    bt_conf <- within(list(), {
	        source(conf.file, local = TRUE)
	    })
	    myargs <- list(
	        cluster = "lsf",
	        registryargs = registryargs,
	        ...
	    )
	    param <- BatchtoolsParam(
	        cluster = "lsf",
	        registryargs = registryargs,
	        ...
	    )
	    if (!is.null(bt_conf$default.resources)) {
	        merged_resources <- bt_conf$default.resources
	        merged_resources[names(param$resources)] <- param$resources
	        param$resources <- merged_resources
	    }
	    bpstart(param)
	    if (!is.null(bt_conf$cluster.functions)) {
	        param$registry$cluster.functions <- bt_conf$cluster.functions
	    }
	    param
	}

  add_decimal_numbers <- function(target_values, value_list) {
	  # Find indices of target values in the value list
	  target_indices <- match(value_list, target_values)
	  
	  # Generate modified values for target values
	  modified_list <- ifelse(
	    !is.na(target_indices),
	    paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
	    value_list
	  )
	  
	  return(modified_list)
	}

duplicate_rows_based_on_count <- function(dataset, column_to_duplicate, count_dataset) {
	    ###this function takes input dataset = dataset you want to modify, column_to_duplicate = the column that determines duplication,   count_dataset = the dataset youre matching. Output of this is dataset modified to have the extra rows to match count_dataset with incremental .1, etc. 
	    overlap <- intersect(dataset[,column_to_duplicate], count_dataset[,column_to_duplicate])

	    #count_values <- table(count_dataset$Geneid_noVersion) #see how many occurences for each value 
		count_values <- table(count_dataset$Geneid_noVersion[count_dataset$Geneid_noVersion %in% overlap])

	    
	    values_not_in_count <- setdiff(dataset$Geneid_noVersion, names(count_values))

	    dataset_1gene <- dataset[match(names(which(count_values == 1)),dataset$Geneid_noVersion),] #subset of 1s

	    dataset_non_count=dataset[match(values_not_in_count,dataset$Geneid_noVersion),]#subset of ids not in T2T

	  #dataset_more <- dataset[names(which(count_values > 1)),]

	  dataset_more <- dataset[match(names(which(count_values > 1)),dataset$Geneid_noVersion),] #subset of 1s
	  assert_that(sum(is.na(dataset_more))==0)
	  #dataset_dup <- data.frame()
	  for (i in seq_len(nrow(dataset_more))) {
	    value <- dataset_more[i, column_to_duplicate] # pick a value
	      occurrence_count <- count_values[value] #count how many times it appears

	      dup <- dataset_more[rep(value, occurrence_count), ]
	      rownames(dup) <- dup$Geneid_noVersion <- paste0(value, ".",seq_len(occurrence_count))
	      if (i == 1) {
	        dataset_dup <- dup
	        }else {
	        dataset_dup <- rbind(dataset_dup, dup)
	        }
	    
	    }
	    rbind(dataset_1gene, dataset_dup, dataset_non_count)
	}


##load necessary files
	genedata_org_hg18=fread("allcount_matrix_2023-08-29.txt", data.table=FALSE) #counts
	rownames(genedata_org_hg18)=genedata_org_hg18$Geneid
	genedata_org_hg18$Geneid=NULL
	dim(genedata_org_hg18)
	genedata_org_hg19=fread("allcount_matrix_2023-08-29.txt", data.table=FALSE) #counts
	rownames(genedata_org_hg19)=genedata_org_hg19$Geneid
	genedata_org_hg19$Geneid=NULL
	dim(genedata_org_hg19)
	genedata_org_v30=fread("allcount_matrix_2023-08-29.txt", data.table=FALSE) #counts
	rownames(genedata_org_v30)=genedata_org_v30$Geneid
	genedata_org_v30$Geneid=NULL
	dim(genedata_org_v30)
	genedata_org_v43=fread("allcount_matrix_2023-08-29.txt", data.table=FALSE) #counts
	rownames(genedata_org_v43)=genedata_org_v43$Geneid
	genedata_org_v43$Geneid=NULL
	dim(genedata_org_v43)
	genedata_org_T2T=fread("allcount_matrix_2023-08-29.txt", data.table=FALSE) #counts
	rownames(genedata_org_T2T)=genedata_org_T2T$Geneid
	genedata_org_T2T$Geneid=NULL
	dim(genedata_org_T2T)

	map <- fread("map_genename_T2T_ensembl_from_gtf.txt") #see AD-DEGs-genome-updates/files/map_genename_T2T_ensembl_from_gtf.txt

	ids_to_keep=readRDS("ids_to_keep.RDS") #see AD-DEGs-genome-updates/files/ids_to_keep.csv

	info_all_base=readRDS("infoall_09.08.23.RDS") #meta info

	main_path="/interaction/"

	all_var <- read.delim("all_var.txt", sep = "\t", header = F)
# 	structure(list(V1 = c("genedata_org_hg18", "genedata_org_hg18", 
# "genedata_org_hg18", "genedata_org_hg18", "genedata_org_hg19", 
# "genedata_org_hg19", "genedata_org_hg19", "genedata_org_v30", 
# "genedata_org_v30", "genedata_org_v43"), V2 = c("hg18", "hg18", 
# "hg18", "hg18", "hg19", "hg19", "hg19", "v30", "v30", "v43"), 
#     V3 = c("genedata_hg18", "genedata_hg18", "genedata_hg18", 
#     "mod_genedata_hg18", "genedata_hg19", "genedata_hg19", "mod_genedata_hg19", 
#     "genedata_v30", "mod_genedata_v30", "mod_genedata_v43"), 
#     V4 = c("genedata_org_hg19", "genedata_org_v30", "genedata_org_v43", 
#     "genedata_org_T2T", "genedata_org_v30", "genedata_org_v43", 
#     "genedata_org_T2T", "genedata_org_v43", "genedata_org_T2T", 
#     "genedata_org_T2T"), V5 = c("hg19", "v30", "v43", "T2T", 
#     "v30", "v43", "T2T", "v43", "T2T", "T2T"), V6 = c("genedata_hg19", 
#     "genedata_v30", "genedata_v43", "genedata_T2T", "genedata_v30", 
#     "genedata_v43", "genedata_T2T", "genedata_v43", "genedata_T2T", 
#     "genedata_T2T"), V7 = c("genedata_org_hg18_hg19", "genedata_org_hg18_v30", 
#     "genedata_org_hg18_v43", "genedata_org_hg18_T2T", "genedata_org_hg19_v30", 
#     "genedata_org_hg19_v43", "genedata_org_hg19_T2T", "genedata_org_v30_v43", 
#     "genedata_org_v30_T2T", "genedata_org_v43_T2T"), V8 = c("hg19_hg18", 
#     "v30_hg18", "v43_hg18", "T2T_hg18", "v30_hg19", "v43_hg19", 
#     "T2T_hg19", "v30_v43", "T2T_v30", "T2T_v43"), V9 = c("info_all_hg18_hg19", 
#     "info_all_hg18_v30", "info_all_hg18_v43", "info_all_hg18_T2T", 
#     "info_all_hg19_v30", "info_all_hg19_v43", "info_all_hg19_T2T", 
#     "info_all_v30_v43", "info_all_v30_T2T", "info_all_v43_T2T"
#     )), class = "data.frame", row.names = c(NA, -10L))

	library(foreach)
	library(doParallel)
	registerDoParallel(cores=nrow(all_var))

	date=""
	path_to_conf_file="/.batchtools.conf.R"

res_final=foreach(i = seq_len(nrow(all_var))) %dopar% {
	if ((all_var[i,3] %in% c("genedata_hg18", "genedata_hg19","genedata_v30", "genedata_v43"))){
		if(length(Sys.glob(paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS")))==0){
		##load first
			genedata_org <- get(all_var[i,4])
			genedata=genedata_org[,6:ncol(genedata_org)]
			##remove .
			genedata$Geneid_noVersion=unlist(lapply(strsplit(as.character(rownames(genedata)),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))
			rownames(genedata) <- genedata$Geneid_noVersion
			genedata_new=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			colnames(genedata_new)=paste(colnames(genedata_new),".", all_var[i,5], sep="")
			genedata_new=as.data.table(genedata_new, keep.rownames=TRUE)
			assign(all_var[i,6],genedata_new)
			print(all_var[i,6])
		##load second 
			genedata_org <- get(all_var[i,1])
			genedata=genedata_org[,6:ncol(genedata_org)]
			genedata$Geneid_noVersion=unlist(lapply(strsplit(as.character(rownames(genedata)),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))
			genedata <- as.data.frame(genedata) 
			rownames(genedata) <- genedata$Geneid_noVersion

			genedata=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			colnames(genedata)=paste(colnames(genedata), ".",all_var[i,2], sep="")
			genedata=as.data.table(genedata, keep.rownames=TRUE)
			assign(all_var[i,3], genedata)
			print(all_var[i,3])
		##merge
			merged_files <- merge(get(all_var[i,3]), get(all_var[i,6]), all = TRUE, by = "rn")
			merged_files<- as.data.frame(merged_files)
			dim(merged_files)
			merged_files[is.na(merged_files)] <- 0
			rownames(merged_files) <- merged_files$rn
			x <- merged_files
			x$rn=NULL
			assign(all_var[i,7],x)
			saveRDS(get(all_var[i,7]), paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS"))
			 }
		else{
			 	x=readRDS(paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS"))
			 	assign(all_var[i,7],x)
			 }
	}
	else{
		##T2T
		if(length(Sys.glob(paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS")))==0){
			genedata_org <- get(all_var[i,4])
			genedata=genedata_org[,6:ncol(genedata_org)]
			genedata$gene_id <- rownames(genedata)
			assert_that(nrow(genedata)==nrow(map))
			genedata <- merge(genedata, map, by.x = "gene_id" , by.y = "gene_id")
			##remove .
			genedata$Geneid_noVersion=unlist(lapply(strsplit(as.character(genedata$projection_parent_gene),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))
			##see which genes are duplicated 
			dupl <- genedata[duplicated(genedata$projection_parent_gene),]

			#unique ids
			uni_ids <- unique(dupl$projection_parent_gene)
			uni_ids <- as.data.frame(uni_ids)
			colnames(uni_ids) <- "ids"
			uni_ids$Geneid_noVersion=unlist(lapply(strsplit(as.character(uni_ids$ids),".",fixed=TRUE),function(x){x[1]}))

			#add decimal 
			values <- genedata$Geneid_noVersion  # Vector of values to modify
			targets <- uni_ids$Geneid_noVersion  # Vector of target values

			modified_values <- add_decimal_numbers(targets, values)
			genedata$modified_values <- unlist(modified_values)

			#make rownames 
			rownames(genedata) <- genedata$modified_values

			genedata_new=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			colnames(genedata_new)=paste(colnames(genedata_new),".", all_var[i,5], sep="")
			genedata_formod=genedata
			genedata_new=as.data.table(genedata_new, keep.rownames=TRUE)
			assign(all_var[i,6],genedata_new)
			print(all_var[i,6])
		#go back to mod ones
			genedata_org <- get(all_var[i,1])
			genedata=genedata_org[,6:ncol(genedata_org)]
			genedata$Geneid_noVersion=unlist(lapply(strsplit(as.character(rownames(genedata)),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))
			genedata <- as.data.frame(genedata) 
			column_to_duplicate <- "Geneid_noVersion"
			duplicated <- duplicate_rows_based_on_count(genedata, column_to_duplicate, genedata_formod)
			rownames(duplicated) <- duplicated$Geneid_noVersion
			#duplicated <- duplicated[, -which(names(duplicated) %in% c("Geneid_noVersion"))] #removes the character column
			genedata=duplicated
			genedata=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			colnames(genedata)=paste(colnames(genedata), ".",all_var[i,2], sep="")
			genedata=as.data.table(genedata, keep.rownames=TRUE)
			assign(all_var[i,3], genedata)
			print(all_var[i,3])
		#merge
			merged_files <- merge(get(all_var[i,3]), get(all_var[i,6]), all = TRUE, by = "rn")
			merged_files<- as.data.frame(merged_files)
			dim(merged_files)
			merged_files[is.na(merged_files)] <- 0
			rownames(merged_files) <- merged_files$rn
			x <- merged_files
			x$rn=NULL
			assign(all_var[i,7],x)
			saveRDS(get(all_var[i,7]), paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS"))
		}else{
		 	x=readRDS(paste0(main_path, all_var[i,8],"/expression/allcount_matrix",sep="_",all_var[i,8], sep="_", date, ".RDS"))
		 	assign(all_var[i,7],x)
		 }

	}
		#covariate
			print(all_var[i,7])
			info_all=info_all_base
			info_all2=info_all_base
			info_all$sample <- paste0(info_all$sample, ".",all_var[i,2], sep="") 
			info_all2$sample <- paste0(info_all2$sample, ".",all_var[i,5], sep="")
			info_all$assembly <-  all_var[i,2]
			info_all2$assembly <-  all_var[i,5]
			info_all3 <- rbind(info_all, info_all2)
			# table(info_all3$assembly)
			assign(all_var[i,9], info_all3)
			#saveRDS(get(all_var[i,9]),paste0(main_path, all_var[i,8],"/covariate/info_all",sep="_",all_var[i,8], sep="_", date, ".RDS"))
			print(all_var[i,9])

		#run dream
			genedata <- get(all_var[i,7])
			isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata) #MAKE SURE AT LEAST 10% OF SAMPLES WITH 1CPM PER GENE
			genesAll <- DGEList(counts=genedata[isexpr,]) #MAKE INTO DGELIST AND REMOVE LOW EXPRESSED GENES
			genesAll <- calcNormFactors(genesAll) #ESTIMATE EFFECTIVE LIBRARY SIZES AND CALCULATE NORMALIZATION FACTORS 
			#dim(genesAll) #checks

			##add covariate info
			info_all3 <- get(all_var[i,9])
			# assert_that(length(setdiff(make.names(as.character(info_all3$sample)), colnames(genesAll)))==0)
			# assert_that(length(setdiff(colnames(genesAll),make.names(as.character(info_all3$sample)))==0))

			# Get the values to match
			values_to_match <- make.names(as.character(info_all3$sample))

			# Identify columns to keep based on matching condition
			columns_to_keep <- colnames(genesAll) %in% values_to_match

			# Subset the data frame to keep only the matching columns
			matched_df <- genesAll[, columns_to_keep]

			# Subset the data frame to remove the non-matching columns
			df_without_nonmatching <- genesAll[, columns_to_keep]

			re_order <- match(make.names(as.character(info_all3$sample)), colnames(df_without_nonmatching))
			colnames(df_without_nonmatching) <- colnames(df_without_nonmatching)[re_order]
			#checks
			assert_that(identical(colnames(df_without_nonmatching),make.names(as.character(info_all3$sample))), msg=paste("assertion 1 on iteration ",i)) #checks

			genesAll <- df_without_nonmatching

			assert_that(length(setdiff(make.names(as.character(info_all3$sample)), colnames(genesAll)))==0, msg=paste("assertion 2 on iteration ",i))
			assert_that(length(setdiff(colnames(genesAll),make.names(as.character(info_all3$sample))))==0, msg=paste("assertion 3 on iteration ",i))

			info_all = info_all3
			info_all$id = info_all$synapse_id
			rownames(info_all)=make.names(as.character(info_all$sample))

			assert_that(identical(colnames(genesAll),rownames(info_all)), msg=paste("assertion 4 on iteration ",i)) #checks
				
			setwd(paste(main_path,all_var[i,8],sep = ""))
			system("mkdir plots")
			pdf("plots/vobj_mean_variance_plot_1cpm_0.1percent.pdf")
			vobj=voom(genesAll,plot=TRUE) # important
			dev.off()

			##formula
			formVoom = formula(paste0("~ assembly + PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)")) 
			# ~assembly + (1 | Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + 
			#     (1 | msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + 
			#     pmi + (1 | race)

			forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

			info_all2=info_all

			##make matrix
			design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formVoom)[2],fixed=TRUE))),info_all2)
			info_all2=info_all2[rownames(design),]
			genedata=genedata[,rownames(design)]
			assert_that(identical(colnames(genedata),make.names(as.character(info_all2$sample))), msg=paste("assertion 5 on iteration ",i))
			 #checks

			isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
			genesAll <- DGEList(counts=genedata[isexpr,])#, group=info$disease) 
			genesAll <- calcNormFactors(genesAll) 
			dim(genesAll) #checks 

			info_all2$id <- as.factor(info_all2$id)
			info_all2$sample <- as.factor(info_all2$sample)
			info_all2$assembly <- as.factor(info_all2$assembly)
			# run voomwithdreamweights
			nonvoom0.5 <- cpm(genesAll,log = TRUE, prior.count = 0.5)
			info_all2$assembly=factor(info_all2$assembly, levels=c("hg18","hg19", "v30", "v43", "T2T"))
			info_all2$assembly=droplevels(info_all2$assembly)

			#more set up 
			#subtrait info
			info_all=info_all2
  
			info_all$CDR_simplified=as.character(info_all$CDR)
			info_all$CDR_simplified[as.numeric(as.character(info_all$CDR))<1]="NL"
			info_all$CDR_simplified[as.numeric(as.character(info_all$CDR))>=1]="AD"
			info_all$CDR_simplified=factor(info_all$CDR_simplified,levels=c("NL","AD"))

			info_all$CERJ_defvsctl=substring(as.character(info_all$CERJ),1,1)
			info_all$CERJ_defvsctl[as.numeric(substring(as.character(info_all$CERJ),1,1))==1]="NL"
			info_all$CERJ_defvsctl[as.numeric(substring(as.character(info_all$CERJ),1,1))==2]="AD"
			info_all$CERJ_defvsctl=factor(info_all$CERJ_defvsctl)

			identical(as.character(info_all$sample),colnames(genesAll))
			info_all2=info_all
			disease_statuses=c("CERJ_defvsctl","CDR_simplified","PlaqueMean")
		for (status in disease_statuses){

			fitDream=NULL
		    info_alltmp=info_all2[which(is.na(info_all2[,status])==F),]
			    resVPtmp=nonvoom0.5[,which(is.na(info_all2[,status])==F)]
			    for(i in 1:ncol(info_alltmp)){
			        if(class(info_alltmp[,i])=="factor"){
			            info_alltmp[,i]=factor(info_alltmp[,i])
			        }
			    }
			# important
			    if(status=="PlaqueMean"){
			         formula_test=formula(paste("~ ", status," * assembly + PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))
			    }else{
			        formula_test=formula(paste("~ 0 + ", status," * assembly + PMI + (1|RACE) + (1|correct_SEX) + RIN + Exonic.Rate + (1|Batch_for_correct)",sep=""))
			    }
			# important
			    design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
			    vobjDreamSubset=nonvoom0.5[,match(rownames(design),colnames(nonvoom0.5))]
			    fullInfo2Subset=info_alltmp[match(rownames(design),make.names(as.character(info_alltmp$sample))),]
			    assert_that(identical(make.names(as.character(fullInfo2Subset$sample)),colnames(vobjDreamSubset)))

			        forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formula_test)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
			        if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
			            fitDream=tryCatch({dream(vobjDreamSubset, formula_test, fullInfo2Subset, BPPARAM = MyBatchtoolsLSFParam(conf.file = "/hpc/users/lunda02/.batchtools.conf.R", workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = "15GB")))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror while running dream on ",as.character(formula_test),": "),e); NULL}) #, ddf="Kenward-Roger") # important
			            if(is.null(fitDream)==FALSE){
			              system(paste0("mkdir DE/",forFileName))
			              saveRDS(fitDream,file=paste0("DE/",forFileName,"/fitDream.RDS"))
			            }
			        }else{
			              fitDream=readRDS(paste0("DE/",forFileName,"/fitDream.RDS"))
			        } 
			       #eBfit=eBayes(fitDream)     # important
			    	coef_regex="PlaqueMean:assembly|NP.14:assembly|AD:assembly|PlaqueMean_defvsctlNL:assembly"
		        	co=grep(coef_regex, colnames(design), value = T) 
		        	assert_that(length(co)==1)
			        res1 = topTable(fitDream, coef=co, number=Inf, sort.by="P")

			        len=length(which(res1$adj.P.Val<0.05))
			    	cat("disease status",status,"DE genes",len,"\n")
			        fwrite(res1,paste("DE/DE_AD_corrected_for_Batch_RINcontinuous_sex_exonicRate_pmi_race_",status,"_test.txt",sep=""),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
	}
}
