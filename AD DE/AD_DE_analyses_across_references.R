##make files needed
	info_all_base=readRDS("info_all.RDS") #meta info
	ids_to_keep=readRDS("voom_normalized_counts_allRin_allpmi_nomislabelled.RDS") ##samples to keep 

	info_all_base <- info_all_base[match(colnames(ids_to_keep),make.names(info_all_base$sample)),]


	info_all_base<- info_all_base[order(info_all_base$sample), ]

	info_all_base$original_sample=info_all_base$sample

	ids_to_keep <- info_all_base$sample

	saveRDS(ids_to_keep,"ids_to_keep.RDS" )
	saveRDS(info_all_base,"infoall_09.13.23.RDS" )

###########################libraries

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

	##function to add decimal - for T2T
	add_decimal_numbers <- function(target_values, value_list) {
	  # Find indices of target values in the value list
	  target_indices <- match(value_list, target_values)
	  
	  # Generate modified values for target values based on how many times it appears
	  modified_list <- ifelse(
	    !is.na(target_indices), #make sure to remove NAs 
	    paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)), #ave produce sequence of #s in order of appearance in valuelist
	    value_list 
	  )
	  
	  return(modified_list)
	}

	duplicate_rows_based_on_count <- function(dataset, column_to_duplicate, count_dataset) {
	    ###this function takes input dataset = dataset you want to modify, column_to_duplicate = the column that determines duplication, count_dataset = the dataset youre matching. Output of this is dataset modified to have the extra rows to match count_dataset with incremental .1, etc. 
	    overlap <- intersect(dataset[,column_to_duplicate], count_dataset[,column_to_duplicate])
	   dataset <- dataset[dataset[,column_to_duplicate] %in% overlap ,]
	   count_dataset <- count_dataset[count_dataset[,column_to_duplicate] %in% overlap ,]
	  
	  count_values <- table(count_dataset[, column_to_duplicate]) #see how many occurences for each value 

	  dataset_1gene <- dataset[names(which(count_values==1)),] #subset of 1s

	  dataset_more <- dataset[names(which(count_values > 1)),]
	  dataset_more <- na.omit(dataset_more) 
	  dataset_dup <- data.frame()
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
	    rbind(dataset_1gene, dataset_dup)
	}

##load raw counts
	genedata_org_hg18=fread("allcount_matrix_2023-04-11.txt", data.table=FALSE) #counts
	rownames(genedata_org_hg18)=genedata_org_hg18$Geneid
	genedata_org_hg18$Geneid=NULL
	dim(genedata_org_hg18)

	genedata_org_hg19=fread("allcount_matrix_2023-04-12.txt", data.table=FALSE) #counts
	rownames(genedata_org_hg19)=genedata_org_hg19$Geneid
	genedata_org_hg19$Geneid=NULL
	dim(genedata_org_hg19)

	genedata_org_v30=fread("allcount_matrix_2023-04-11.txt", data.table=FALSE) #counts
	rownames(genedata_org_v30)=genedata_org_v30$Geneid
	genedata_org_v30$Geneid=NULL
	dim(genedata_org_v30)

	genedata_org_v43=fread("allcount_matrix_2023-04-11.txt", data.table=FALSE) #counts
	rownames(genedata_org_v43)=genedata_org_v43$Geneid
	genedata_org_v43$Geneid=NULL
	dim(genedata_org_v43)

	genedata_org_T2T=fread("allcount_matrix_2023-07-03.txt", data.table=FALSE) #counts
	rownames(genedata_org_T2T)=genedata_org_T2T$Geneid
	genedata_org_T2T$Geneid=NULL
	dim(genedata_org_T2T)

	map <- fread("map_genename_T2T_ensembl_from_gtf.txt") #see AD-DEGs-genome-updates/files/map_genename_T2T_ensembl_from_gtf.txt

	all_var <- 	read.delim("/sc/arion/projects/mscic1/results/anina/fun_project_4.23/run_dge_AD.txt", header = FALSE)

# 	structure(list(V1 = c("genedata_org_hg18", "genedata_org_hg19", 
# "genedata_org_v30", "genedata_org_v43", "genedata_org_T2T"), 
#     V2 = c("info_all_hg18", "info_all_hg19", "info_all_v30", 
#     "info_all_v43", "info_all_T2T"), V3 = c("RosmapHG18", "RosmapHG19", 
#     "Rosmapv30", "Rosmapv43", "RosmapT2T"), V4 = c("hg18", "hg19", 
#     "v30", "v43", "T2T")), class = "data.frame", row.names = c(NA, 
# -5L))


	ids_to_keep=readRDS("ids_to_keep.RDS") #see AD-DEGs-genome-updates/files/ids_to_keep.csv

	info_all_base=readRDS("infoall_09.13.23.RDS") #meta info

	main_path= #where to save

	library(foreach)
	library(doParallel)
	registerDoParallel(cores=nrow(all_var))

	date= #date to save
	path_to_conf_file=".batchtools.conf.R"

###run all
	res_final=foreach(i = seq_len(nrow(all_var))) %dopar% {
	if ((all_var[i,1] %in% c("genedata_org_hg18", "genedata_org_hg19","genedata_org_v30", "genedata_org_v43"))){
		if(length(Sys.glob(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS")))==0){
		##load
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
			# dim(genedata)
			rownames(genedata)=genedata$Geneid_noVersion
			genedata=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			assign(all_var[i,1],genedata)
			print(all_var[i,1])
			saveRDS(get(all_var[i,1]),paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
		}
		else{
			 	x=readRDS(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
			 	assign(all_var[i,1],x)
			 }
			}
	else{
		##T2T
		if(length(Sys.glob(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS")))==0){
			genedata_org <- get(all_var[i,1])
			genedata=genedata_org[,6:ncol(genedata_org)]
			genedata$gene_id <- rownames(genedata)
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
			uni_ids$Geneid_noVersion=unlist(lapply(strsplit(as.character(uni_ids$ids),".",fixed=TRUE),function(x){
				if(sum(grepl("_", x))==FALSE){
					x[1]
				}else{
					paste(x[1], paste(unlist(strsplit(as.character(x[2]),"_"))[2:3],collapse = "_"), sep="_")
				}
				}))

			#add decimal 
			values <- genedata$Geneid_noVersion  # Vector of target values
			targets <- uni_ids$Geneid_noVersion  # Vector of values to modify

			modified_values <- add_decimal_numbers(targets, values)
			genedata$modified_values <- unlist(modified_values)

			#make rownames 
			rownames(genedata) <- genedata$modified_values
			genedata_new=genedata[, intersect(colnames(genedata), make.names(as.character(ids_to_keep)))]
			assign(all_var[i,1],genedata_new)
			print(all_var[i,1])
			saveRDS(get(all_var[i,1]),paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
		}else{
		 	x=readRDS(paste0(main_path, all_var[i,3],"/expression/allcount_matrix",sep="_",all_var[i,4], sep="_", date, ".RDS"))
		 	assign(all_var[i,1],x)
		 }}

		#covariate
			print(all_var[i,2])
			info_all=info_all_base
			assign(all_var[i,2], info_all)
			saveRDS(get(all_var[i,2]),paste0(main_path, all_var[i,3],"/covariate/info_all",sep="_",all_var[i,4], sep="_", date, ".RDS"))
			print(all_var[i,2])

		#run dream
			genedata <- get(all_var[i,1])
			isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata) #MAKE SURE AT LEAST 10% OF SAMPLES WITH 1CPM PER GENE
			genesAll <- DGEList(counts=genedata[isexpr,]) #MAKE INTO DGELIST AND REMOVE LOW EXPRESSED GENES
			genesAll <- calcNormFactors(genesAll) #ESTIMATE EFFECTIVE LIBRARY SIZES AND CALCULATE NORMALIZATION FACTORS 
			#dim(genesAll) #checks

			##add covariate info
			info_all3 <- get(all_var[i,2])
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
			df_without_nonmatching <- df_without_nonmatching[,re_order]
			#checks
			assert_that(identical(colnames(df_without_nonmatching),make.names(as.character(info_all3$sample))), msg=paste("assertion 1 on iteration ",i)) #checks

			genesAll <- df_without_nonmatching

			assert_that(length(setdiff(make.names(as.character(info_all3$sample)), colnames(genesAll)))==0, msg=paste("assertion 2 on iteration ",i))
			assert_that(length(setdiff(colnames(genesAll),make.names(as.character(info_all3$sample))))==0, msg=paste("assertion 3 on iteration ",i))

			info_all = info_all3
			info_all$id = info_all$original_sample
			rownames(info_all)=make.names(as.character(info_all$sample))

			assert_that(identical(colnames(genesAll),rownames(info_all)), msg=paste("assertion 4 on iteration ",i)) #checks
				
			setwd(paste(main_path,all_var[i,3],sep = ""))
			system("mkdir plots")
			pdf("plots/vobj_mean_variance_plot_1cpm_0.1percent.pdf")
			vobj=voom(genesAll,plot=TRUE) # important
			dev.off()

			##formula
			formVoom = formula(paste0("~ (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)")) 
			#~(1 | Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + 
			#    (1 | msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + 
			#    pmi + (1 | race)

			forFileNameVoom=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formVoom)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))

			info_all2=info_all

			##make matrix
			design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formVoom)[2],fixed=TRUE))),info_all2)
			info_all2=info_all2[rownames(design),]
			genedata=genedata[,rownames(design)]
			assert_that(identical(colnames(genedata),make.names(as.character(info_all2$sample)))) 

			isexpr = rowSums(cpm(genedata)>=1) >= 0.1*ncol(genedata)
			genesAll <- DGEList(counts=genedata[isexpr,])#, group=info$disease) 
			genesAll <- calcNormFactors(genesAll) 
			dim(genesAll) #checks 
			assert_that(identical(colnames(genesAll),make.names(as.character(info_all2$sample)))) #checks 

			info_all2$sample <- as.factor(info_all2$sample)

			# run voomwithdreamweights
			if(length(which(paste0(forFileNameVoom)==list.files("expression/")))==0){
			  system(paste0("mkdir -p expression/",forFileNameVoom))
			  # vobjDream = voomWithDreamWeights( genesAll, formVoom, info_all2, save.plot=TRUE, BPPARAM = MulticoreParam(nprocs))
			  vobjDream = tryCatch({voomWithDreamWeights( genesAll, formVoom, info_all2, save.plot=TRUE, BPPARAM = MyBatchtoolsLSFParam(conf.file = path_to_conf_file, workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = "15GB")))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror during voomWithDreamWeights for iteration",as.character(i),": "),e); NULL})
			  saveRDS(vobjDream,file=paste0("expression/",forFileNameVoom,"/voomWithDreamWeights.RDS"))
			}else{
			  vobjDream=readRDS(paste0("expression/",forFileNameVoom,"/voomWithDreamWeights.RDS"))
			}

			#subtrait info
			info_all2$braaksc_simplified="NL"
			info_all2$braaksc_simplified[as.numeric(as.character(info_all2$braaksc))>=3]="AD"
			info_all2$braaksc_simplified=factor(info_all2$braaksc_simplified,levels=c("NL","AD"))

			info_all2$ceradsc_simplified="NL"
			info_all2$ceradsc_simplified[as.numeric(as.character(info_all2$ceradsc))<4]="AD"
			info_all2$ceradsc_simplified=factor(info_all2$ceradsc_simplified,levels=c("NL","AD"))

			info_all2$cogdx_simplified="NL"
			info_all2$cogdx_simplified[as.numeric(as.character(info_all2$cogdx))>=4]="AD"
			info_all2$cogdx_simplified[is.na(info_all2$cogdx)]=NA
			info_all2$cogdx_simplified=factor(info_all2$cogdx_simplified,levels=c("NL","AD"))

			info_all2$cogdx_defvsctl="MCI"
			info_all2$cogdx_defvsctl[as.numeric(as.character(info_all2$cogdx))==1]="NL"
			info_all2$cogdx_defvsctl[as.numeric(as.character(info_all2$cogdx))==4 | as.numeric(as.character(info_all2$cogdx))==5]="AD"
			info_all2$cogdx_defvsctl[is.na(info_all2$cogdx)]=NA
			info_all2$cogdx_defvsctl=factor(info_all2$cogdx_defvsctl,levels=c("NL","AD","MCI"))

			disease_statuses=c("cts_mmse30_lv","braaksc_simplified","ceradsc_simplified","cogdx_simplified","ceradsc","cogdx_defvsctl")

 			for (status in disease_statuses){

			fitDream=NULL
			if(status!="ceradsc_defvsctl"){
		        info_alltmp=info_all2[which(is.na(info_all2[,status])==F),]
		        resVPtmp=vobjDream[,which(is.na(info_all2[,status])==F)]
		    }else{
		        info_alltmp=info_all2[which(is.na(info_all2[,"ceradsc"])==F),]
		        resVPtmp=vobjDream[,which(is.na(info_all2[,"ceradsc"])==F)]
		    }
		    for(i in 1:ncol(info_alltmp)){
		        if(class(info_alltmp[,i])=="factor"){
		            info_alltmp[,i]=factor(info_alltmp[,i])
		        }
		    }

		    if(status=="cts_mmse30_lv"){
		        formula_test=formula(paste("~ ",status," + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))
		    }else if(status!="ceradsc_defvsctl"){
		        formula_test=formula(paste("~ 0 + ",status," + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))
		    }else{
		        formula_test=formula(paste("~ 0 + ceradsc + (1|Batch) + RINcontinuous + RnaSeqMetrics__MEDIAN_5PRIME_TO_3PRIME_BIAS + (1|msex) + AlignmentSummaryMetrics__STRAND_BALANCE + RnaSeqMetrics__PCT_INTRONIC_BASES + pmi + (1|race)",sep=""))
		    }

			# important
			    design=model.matrix(formula(paste0("~",gsub("1 |","",as.character(formula_test)[2],fixed=TRUE))),info_alltmp)
		    	vobjDreamSubset=vobjDream[,match(rownames(design),colnames(vobjDream))]
		    	fullInfo2Subset=info_alltmp[match(rownames(design),make.names(as.character(info_alltmp$sample))),]
			    assert_that(all(make.names(as.character(fullInfo2Subset$sample))==colnames(vobjDreamSubset)), msg=paste("assertion 6 on iteration ",i))

			    if(status=="cts_mmse30_lv"){
		        forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formula_test)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
		        if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
		            fitDream=tryCatch({dream(vobjDreamSubset, formula_test, fullInfo2Subset, BPPARAM = MyBatchtoolsLSFParam(conf.file = path_to_conf_file, workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = "15GB")))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror while running dream on ",as.character(formula_test),": "),e); NULL}) #, ddf="Kenward-Roger") # important
		            if(is.null(fitDream)==FALSE){
		              system(paste0("mkdir DE/",forFileName))
		              saveRDS(fitDream,file=paste0("DE/",forFileName,"/fitDream.RDS"))
		            }
		        }else{
		              fitDream=readRDS(paste0("DE/",forFileName,"/fitDream.RDS"))
		        } 
		        eBfit=eBayes(fitDream)     # important
		        res1 = topTable(eBfit, coef="cts_mmse30_lv", number=Inf, sort.by="P")

		        len=length(which(res1$adj.P.Val<0.05))
		    	cat("assembly for iteration",status,"DE genes",len,"\n")

		    	fwrite(res1,paste("DE/DE_AD_corrected_for_Batch_RINcontinuous_MEDIAN_5PRIME_TO_3PRIME_BIAS_msex_STRAND_BALANCE_PCT_INTRONIC_BASES_pmi_race_",status,"_test.txt",sep=""),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
		    }else{
		 # important
		         if(status=="braaksc_simplified"){
		            L = getContrast( vobjDreamSubset, formula_test, fullInfo2Subset, c("braaksc_simplifiedAD", "braaksc_simplifiedNL"))
		        }else if(status=="ceradsc_simplified"){
		            L = getContrast( vobjDreamSubset, formula_test, fullInfo2Subset, c("ceradsc_simplifiedAD", "ceradsc_simplifiedNL"))
		        }else if(status=="cogdx_simplified"){
		            L = getContrast( vobjDreamSubset, formula_test, fullInfo2Subset, c("cogdx_simplifiedAD", "cogdx_simplifiedNL"))
		        }else if(status=="ceradsc"){
		            L = getContrast( vobjDreamSubset, formula_test, fullInfo2Subset, c(paste0("ceradsc","1"), paste0("ceradsc","4")))
		        }else if(status=="cogdx_defvsctl"){
		            L = getContrast( vobjDreamSubset, formula_test, fullInfo2Subset, c("cogdx_defvsctlAD", "cogdx_defvsctlNL"))
		        }
		 # important

		        forFileName=gsub("Corrected_","",gsub("PERCENT|Percent","PCT",gsub("PRIME","P",gsub("_scaled|RNASEQ_|RNA_|RnaSeqMetrics__|AlignmentSummaryMetrics__","",make.names(gsub("1|","",gsub(")","",gsub("(","",gsub(" ","",gsub(" + ","_",as.character(formula_test)[2],fixed=TRUE)),fixed=TRUE),fixed=TRUE),fixed=TRUE))))))
		        if(length(which(paste0(forFileName,"/fitDream.RDS")==list.files("DE/",recursive=TRUE)))==0){
		            fitDream=tryCatch({dream(vobjDreamSubset, formula_test, fullInfo2Subset, L, BPPARAM = MyBatchtoolsLSFParam(conf.file = "/hpc/users/lunda02/.batchtools.conf.R", workers = 100, resources = list(walltime = 30 * 60, threads_per_job = 1, memory = "15GB")))}, error = function(e){message(paste0("\n\n\n\n\n\n\n\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\n################\nerror while running dream on ",as.character(formula_test),": "),e); NULL}) #, ddf="Kenward-Roger")
		            if(is.null(fitDream)==FALSE){
		              system(paste0("mkdir DE/",forFileName))
		              saveRDS(fitDream,file=paste0("DE/",forFileName,"/fitDream.RDS"))
		            }
		        }else{
		              fitDream=readRDS(paste0("DE/",forFileName,"/fitDream.RDS"))
		        }
		        eBfit = eBayes(fitDream)
		        res1 = topTable(eBfit, coef="L1", number=Inf, sort.by="P")
		    len=length(which(res1$adj.P.Val<0.05))
		    cat("assembly for iteration",status,"DE genes",len,"\n")

		    fwrite(res1,paste("DE/DE_AD_corrected_for_Batch_RINcontinuous_MEDIAN_5PRIME_TO_3PRIME_BIAS_msex_STRAND_BALANCE_PCT_INTRONIC_BASES_pmi_race_",status,"_test.txt",sep=""),col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
		}
		cat(all_var[i,3])
	}
}

			    


