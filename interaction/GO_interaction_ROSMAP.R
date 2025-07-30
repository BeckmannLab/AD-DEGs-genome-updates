library(ggplot2)
options(stringsAsFactors = FALSE)
library("R.matlab")
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)
library(data.table)

library("AnnotationDbi")
library("GO.db")
library("plyr")
library("org.Mm.eg.db")



annotate_GOterms_DE <- function(mel,sample_name,ensGene=TRUE,revigo=FALSE,GOFigure=FALSE,rrvgo=TRUE,gene.map=gene.map, useFDR=FALSE){
    # Performs GO annotation of the genes of interest
    # inputs are mel, a 2 column data.frame with number of rows corresponding
    # to your sample size, where the first column is gene name and the second
    # column is a binary vector of 0s and 1s, with 0s for genes in background
    # (not significantly differentially expressed genes for example) and 1s
    # for genes of interest (significantly differentially expressed genes for
    # example); sample_name is a string containing the name of the file that 
    # will be saved with the output of the analysis; and date, a string of the
    # date of the experiment that will also be used in the name of the output
    # file. The function has no output, output is directly saved as a table.
    #gofigure open source while revigo is not 
    options(stringsAsFactors = FALSE)
    library("R.matlab")
    library(goseq)
    library(topGO)
    library(org.Hs.eg.db)
    library(Rgraphviz)
    library(rrvgo)
    library(Matrix.utils)

    # load in the data set
    #
    #######################################################################################
    #
    # now carry out GO analysis
    #
    # read in analysis results above if needed
    #

    #pwf = nullp(a, 'hg19', 'geneSymbol')

    #rename and output of file - don't want to override mel object 
    net = mel
    output.annotation.file = paste(sample_name,"_module_enrichments",".txt",sep="")
    
    #make directories 
    if(revigo==TRUE){
        output.annotation.revigo = paste(sample_name,"_","revigo",sep="")
        if(dir.exists(output.annotation.revigo)){
          system(paste("rm -rf ",output.annotation.revigo,sep=""))  #remove folder if folder already exists 
        }
        dir.create(output.annotation.revigo,showWarnings=FALSE)
    }
    if(GOFigure==TRUE){
        output.annotation.GOFigure = paste(sample_name,"_","GOFigure",sep="")
        if(dir.exists(output.annotation.GOFigure)){
          system(paste("rm -rf ",output.annotation.GOFigure,sep=""))
        }
        dir.create(output.annotation.GOFigure,showWarnings=FALSE)
    }
    if(rrvgo==TRUE){
        output.annotation.rrvgo = paste(sample_name,"_","rrvgo",sep="")
        if(dir.exists(output.annotation.rrvgo)){
          system(paste("rm -rf ",output.annotation.rrvgo,sep=""))
        }
        dir.create(output.annotation.rrvgo,showWarnings=FALSE)
    }

    #get the GO IDS - gene set. which gene set belong to 
    #commented when running in a parallel loop and do it before you run if doing in parallel 
    if (is.null(gene.map)){
    if(ensGene==TRUE){
        gene.map = getgo(net[,1],'hg19','ensGene')       # extracting go IDs from gene names
    }else{
        gene.map = getgo(net[,1],'hg19','geneSymbol')       # extracting go IDs from gene names
    }
}

    gene.map = gene.map[!is.na(names(gene.map))]        # removing non existing entries
    number_of_GO = length(unique(unlist(gene.map))) #21722
    
    #safty to make sure something not loaded in
    if (exists("res")){
       remove("res")
    }

    counter = 0
    counter = counter + 1
    flush.console()
 
    x = net[,1]         # extracting gene names
    a = rep(0,length(x))    # creating vector of length module names filled with 0s and name of row gene name
    names(a) = x
    a[net[,2] == 1] = 1     # changing indices in this vector at module positions to 1, vector used in definition of GO element

    ##BP = biological processes 
    ips_BP = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("BP"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)
    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_BP, test.stat)   # doing Fisher statistics on GO graph
    res.final_BP = GenTable(ips_BP, classic = res.fisher, topNodes=length(ips_BP@graph@nodes))   # selecting top 200 of res.fisher

    #mf = molecular functinos 
    ips_MF = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("MF"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)


    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_MF, test.stat)   # doing Fisher statistics on GO graph
    res.final_MF = GenTable(ips_MF, classic = res.fisher, topNodes=length(ips_MF@graph@nodes))   # selecting top 200 of res.fisher

    #cc = cellular components 
    ips_CC = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
                     ontology = c("CC"), 
                     allGenes = as.factor(a), 
                     geneSel = names(a[a==1]),
                     nodeSize = 10,
                     annot = annFUN.gene2GO,
                     gene2GO = gene.map)
    test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    res.fisher = getSigGroups(ips_CC, test.stat)   # doing Fisher statistics on GO graph
    res.final_CC = GenTable(ips_CC, classic = res.fisher, topNodes=length(ips_CC@graph@nodes))   # selecting top 200 of res.fisher

    # ips = new("topGOdata", description = "Enrichment in mono cells", # creating GO element and graph
    #                  ontology = c("BP", "MF", "CC", "KEGG"), 
    #                  allGenes = as.factor(a), 
    #                  geneSel = names(a[a==1]),
    #                  nodeSize = 10,
    #                  annot = annFUN.gene2GO,
    #                  gene2GO = gene.map)

    # test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test") # creating element of class classic count
    # res.fisher = getSigGroups(ips, test.stat)   # doing Fisher statistics on GO graph
    # res.final = GenTable(ips, classic = res.fisher, topNodes=length(ips@graph@nodes))   # selecting top 200 of res.fisher

    # gene.list=get_genes_in_signif_pathways_GO_ensembl(res.final[,"GO.ID"],gene.map)

    res=cbind(rbind(res.final_BP,res.final_MF,res.final_CC))


    # if (exists("res")){
    #     res = rbind(res, cbind(res.final))
    # }else{
    #     res = cbind(res.final) # creating res table with 1st column module color and rest res.final
    # }

    #BH = benjamini-hochberg procedure - decreases FDR. avoid type 1 errors


    res.mod = cbind(res, res[,"Significant"]/res[,"Expected"],p.adjust(res[,"classic"],method="BH",n=number_of_GO))  # adding a column of ratio of significant over expected (fold enrichment)
    names(res.mod)[c(7,8)] = c("fold_enrichment","BH")
    res.mod = res.mod[,c( "GO.ID", "Term", "Annotated", "Significant", "Expected", "fold_enrichment", "classic","BH")]

    write.table(res.mod, output.annotation.file, sep="\t", quote=FALSE, row.names=FALSE)

    if (useFDR==TRUE){var="BH"}else{var="classic"}
    
    if(rrvgo==TRUE){
    res.mod[is.na(res.mod[,var]),var] = 10^-300 
    tmp = res.mod[res.mod[,var] <=0.05,c('GO.ID',var)]
    go_analysis <- tmp
    simMatrix_BP <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
    simMatrix_MF <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="MF",
                                method="Rel")
    simMatrix_CC <- calculateSimMatrix(go_analysis$GO.ID,
                                orgdb="org.Hs.eg.db",
                                ont="CC",
                                method="Rel")

    combo_simMat <- rBind.fill(simMatrix_MF, simMatrix_CC)
    combo_all <- rBind.fill(combo_simMat, simMatrix_BP)
    go_analysis[,var] <- as.numeric(go_analysis[,var])
    scores <- setNames(-log10(go_analysis[,var]), go_analysis$GO.ID)
    reducedTerms <- reduceSimMatrix(combo_all,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

    pdf(paste(output.annotation.rrvgo, "treemap_output.annotation_DE_AD.pdf", sep=""))
    scatterPlot(combo_all,reducedTerms)
    wordcloudPlot(reducedTerms,min.freq=1, colors="black" )
    heatmapPlot(combo_all,reducedTerms, annotateParent=TRUE, annotationLabel="parentTerm", fontsize=6 )
    treemapPlot(reducedTerms)
    dev.off()
}
                 
    
    
    if(revigo==TRUE){
        res.mod[is.na(res.mod[,"BH"]),"BH"] = 0 

        tmp = res.mod[res.mod[,"BH"] <=1,c('GO.ID','BH')]
        if(nrow(tmp)>0){
            FN=paste(output.annotation.revigo,sep="")
            write.table(res.mod[res.mod[,"BH"] <=1,c('GO.ID','BH')], sep="\t",file=paste(FN,"/GOterms",sep=""),quote=F,col.names=F,row.names=F)
        }   
    # for i in ~/minerva/coexpression_final/revigo/GOterms*; do echo $i;java -jar RevigoStandalone.jar $i --cutoff=0.7 --stdout >${i}_0.7;java -jar RevigoStandalone.jar $i --cutoff=0.1 --stdout >${i}_0.1; Rscript ~/minerva/scripts/Plot_ReviGO_Modules.R $i;done
        folder=getwd()
        system(paste("module unload java; module load java/1.8.0_66; cd ~/source/RevigoStandalone_2015-02-17_beta/; ","for i in ",folder,"/",output.annotation.revigo,"/GOterms*; do java -Xmx4000m -jar RevigoStandalone.jar $i --cutoff=0.7 --stdout >${i}_0.7; java -Xmx4000m -jar RevigoStandalone.jar $i --cutoff=0.1 --stdout >${i}_0.1; done",sep=""))
        if(length(which(res.mod[,"BH"]<1))>1){
            name=paste(output.annotation.revigo,"/GOterms",sep="")
            Plot_ReviGO_Modules(name)
        }
        unlink(paste(folder,"/",output.annotation.revigo,"/*_0.1",sep=""))
        unlink(paste(folder,"/",output.annotation.revigo,"/*_0.7",sep=""))
    }

    if (useFDR==TRUE){var="BH"}else{var="classic"}
    if(GOFigure==TRUE){
        res.mod[is.na(res.mod[,var]),var] = 10^-300 
        tmp = res.mod[res.mod[,var] <=0.05,c('GO.ID',var)]
        if(nrow(tmp)>0){
            FN=paste(output.annotation.GOFigure,sep="")
            write.table(res.mod[res.mod[,var] <=0.05,c('GO.ID',var)], sep="\t",file=paste(FN,"/GOterms",sep=""),quote=F,col.names=F,row.names=F)
          folder=getwd()
          system(paste0("ml python/3.7.3; python /GOFigure/GO-Figure/gofigure.py -i", folder,"/",output.annotation.GOFigure,"/GOterms -o ",folder,"/",output.annotation.GOFigure,"/ -w GOFigure_0.1 -si 0.1; python /sc/arion/projects/mscic1/data/GOFigure/GO-Figure/gofigure.py -i", folder,"/",output.annotation.GOFigure,"/GOterms -o ",folder,"/",output.annotation.GOFigure,"/ -w GOFigure_0.7 -si 0.7;"))
        }
        # python /GO-Figure/gofigure.py -i $folder/GOterms -o $folder/ -w GOFigure_0.1 -si 0.1
        # .1 AND .7 smaller and bigger group 
        if(length(which(res.mod[,var]<0.05))>1){
            name=paste("GOFigure",sep="")
            Plot_GOFigure_Modules(name,output.annotation.GOFigure)
            unlink(paste(folder,"/",output.annotation.GOFigure,"/GOFigure_0.1",sep=""))
            unlink(paste(folder,"/",output.annotation.GOFigure,"/GOFigure_0.7",sep=""))
        }
    }
}


Plot_GOFigure_Modules <- function(name,output.annotation.GOFigure){
  options(stringsAsFactors=F)
  library(treemap)
  files_toRead_0.1=Sys.glob(paste0(output.annotation.GOFigure,"/*_full_table_GOFigure_0.1.tsv"))
  files_toRead_0.7=Sys.glob(paste0(output.annotation.GOFigure,"/*_full_table_GOFigure_0.7.tsv"))

  gofigure.1=c()
  for(i in 1:length(files_toRead_0.1)){
    gofiguretmp.1=fread(files_toRead_0.1[i],data.table=FALSE)
    gofigure.1=rbind(gofigure.1,gofiguretmp.1)
  }

  gofigure.7=c()
  for(i in 1:length(files_toRead_0.7)){
    gofiguretmp.7=fread(files_toRead_0.7[i],data.table=FALSE)
    gofigure.7=rbind(gofigure.7,gofiguretmp.7)
  }

  # gofigure1.1=fread(paste0(output.annotation.GOFigure,"/biological_process_full_table_",name,"_0.1.tsv"),data.table=FALSE)
  # gofigure2.1=fread(paste0(output.annotation.GOFigure,"/cellular_component_full_table_",name,"_0.1.tsv"),data.table=FALSE)
  # gofigure3.1=fread(paste0(output.annotation.GOFigure,"/molecular_function_full_table_",name,"_0.1.tsv"),data.table=FALSE)
  # gofigure1.7=fread(paste0(output.annotation.GOFigure,"/biological_process_full_table_",name,"_0.7.tsv"),data.table=FALSE)
  # gofigure2.7=fread(paste0(output.annotation.GOFigure,"/cellular_component_full_table_",name,"_0.7.tsv"),data.table=FALSE)
  # gofigure3.7=fread(paste0(output.annotation.GOFigure,"/molecular_function_full_table_",name,"_0.7.tsv"),data.table=FALSE)
  # gofigure.1=rbind(gofigure1.1,gofigure2.1,gofigure3.1)
  # gofigure.7=rbind(gofigure1.7,gofigure2.7,gofigure3.7)

  # name_tmp = strsplit(name,"_",fixed=T)[[1]]
  # name = paste(name_tmp[length(name_tmp)-1],"-",name_tmp[length(name_tmp)])
  # by default, outputs to a PDF file
  #pdf( file=paste0(output.annotation.GOFigure,"/",name,"_GOFigure_TreeMap.pdf",sep=""), width=16, height=9 ) # width and height are in inches
pdf(paste(main_path, "GOFigure_TreeMap_DE_AD.pdf", sep=""),width=16, height=9)


  l_c = gofigure.1
  s_c = gofigure.7

  s_c$eliminated=abs(as.numeric(is.element(s_c$`Cluster member`,s_c$`Cluster representative`))-1)
  s_c$real_representative = sapply(s_c$`Cluster member`,function(x){
      l_c$`Cluster representative`[l_c$`Cluster member` == x]
      })
  clusters = s_c[s_c$eliminated!=1,]
  # clusters$cluster_name = sapply(clusters$real_representative,function(X){
  #     clusters$`Cluster member description`[clusters$`Cluster member`==X]
  #     })
  clusters$cluster_name = sapply(clusters$real_representative,function(X){ 
    if(sum(clusters$`Cluster member`==X)>0){
      clusters$`Cluster member description`[clusters$`Cluster member`==X]
      }else{
        ""
      }
      })

  clusters$log10_p.value <- abs(-log10(as.numeric( as.character(clusters$`Member P-value`) )));
  clusters$frequency <- as.numeric( as.character(clusters$`Member frequency`) );
    # check the tmPlot command documentation for all possible parameters - there are a lot more

  treemap(
      clusters,
    index = c("cluster_name","Cluster member description"),
    vSize = "log10_p.value",
    type = "categorical",
    vColor = "cluster_name",
    title = paste(name,"REVIGO Gene Ontology treemap"),
    inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
                                                         # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
      position.legend = "none"
  )

  dev.off()
  }

Plot_ReviGO_Modules <- function(name){
    options(stringsAsFactors=F)
    library(treemap)

    large_clusters = paste(name,"_0.1",sep="") 
    small_clusters = paste(name,"_0.7",sep="")
    outputFileName = paste(name,"_ReviGO.pdf",sep="")
    # name_tmp = strsplit(name,"_",fixed=T)[[1]]
    # name = paste(name_tmp[length(name_tmp)-1],"-",name_tmp[length(name_tmp)])
    # by default, outputs to a PDF file
    pdf(paste(main_path, "ReviGO_DE_AD.pdf"),width=16, height=9 ) # width and height are in inches

    # l_c = read.delim(large_clusters,skip=2,header=T)
    # s_c = read.delim(small_clusters,skip=2, header=T)
    l_c = gofigure.1
    s_c = gofigure.7

    s_c$real_representative = sapply(s_c$term_ID,function(x){
        l_c$representative[l_c$term_ID == x]
        })
    clusters = s_c[s_c$eliminated!=1,]
    clusters$cluster_name = sapply(clusters$real_representative,function(X){
        clusters$description[clusters$term_ID==X]
        })
    # check the tmPlot command documentation for all possible parameters - there are a lot more


    clusters$log10_p.value <- abs(as.numeric( as.character(clusters$log10_p.value) ));
    clusters$frequency <- as.numeric( as.character(clusters$frequency) );
    clusters$uniqueness <- as.numeric( as.character(clusters$uniqueness) );
    clusters$dispensability <- as.numeric( as.character(clusters$dispensability) );

    treemap(
        clusters,
        index = c("cluster_name","description"),
        vSize = "log10_p.value",
        type = "categorical",
        vColor = "cluster_name",
        title = paste(name,"REVIGO Gene Ontology treemap"),
        inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
        lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
        bg.labels = "#CCCCCCAA",     # define background color of group labels
                                                           # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
        position.legend = "none"
    )

    dev.off()
    }

##format
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


df_all_rosmap = readRDS("rosmap_ad_interaction_results_matrix.RDS") #see AD-DEGs-genome-updates/misc/rosmap_ad_interaction_results_matrix.R
common = df_all_rosmap[which(!is.na(df_all_rosmap$DE_AD_sign_same)),]

traits = c("braaksc_simplified", "cogdx_simplified", "ceradsc_defvsctl","cts_mmse30" )
pairs = c("hg19_hg18","v30_hg18","v30_hg18", "v43_hg18","T2T_hg18","v30_hg19", "v43_hg19", "T2T_hg19", "v43_v30", "T2T_v30", "T2T_v43")
map <- fread("/gene_ids_ensembl2symbol_fromHUGO_10JUN2020.tsv")[,.(symbol=`Approved symbol`, gene=`Ensembl gene ID`)] #see AD-DEGs-genome-updates/files/gene_ids_ensembl2symbol_fromHUGO_10JUN2020.csv

all = expand.grid(traits, pairs)

main_path="/interaction_allPathways/"

setwd(main_path)
####setup
##make the files that are input_msbb
for (i in 1:nrow(all)){
    df_sig2 = common[which(common$trait == all[i,1]),]
    df_sig3 = df_sig2[which(df_sig2$assembly_comparison == all[i,2]),]
    for_gsea_df = format_for_gsea(data.frame(gene_symbol = df_sig3$gene_symbol, sig = df_sig3$interaction_status))
    df = for_gsea_df
    df$DE <- ifelse(df$DE == "Signif", 1, 0)
    de2 = df[,c("gene_ID", "DE")]
    gene.map = getgo(unlist(lapply(strsplit(unique(de2$gene_ID),".",fixed=TRUE),function(x){x[1]})),'hg19','geneSymbol') 
    try(annotate_GOterms_DE(mel=de2,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
}



