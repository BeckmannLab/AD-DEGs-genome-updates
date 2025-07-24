rm(list=ls())


library(ggplot2)
options(stringsAsFactors = FALSE)
library("R.matlab")
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)
library(data.table)



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

    res=cbind(rbind(res.final_BP,res.final_MF,res.final_CC))

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
          system(paste0("ml python/3.7.3; python /gofigure.py -i", folder,"/",output.annotation.GOFigure,"/GOterms -o ",folder,"/",output.annotation.GOFigure,"/ -w GOFigure_0.1 -si 0.1; python /gofigure.py -i", folder,"/",output.annotation.GOFigure,"/GOterms -o ",folder,"/",output.annotation.GOFigure,"/ -w GOFigure_0.7 -si 0.7;"))
        }
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


##run
all_results=readRDS("all_resultsDE_1.08.24.RDS") #AD-DEGs-genome-updates/misc/all_resultsDE_1.08.24.RDS

main_path="/allPathways/"
system(paste("mkdir",main_path))
setwd(main_path)

phenos=c("braaksc_simplified","ceradsc_defvsctl","cogdx_simplified","cts_mmse30")
assemblies=c("hg18","hg19","v30","v43","T2T")

assembliesComb = t(combn(assemblies,2))
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[2],x[1])}))

grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")


for(i in 1:nrow(grid)){
    outCol=paste0(grid[i,"assembly2"],grid[i,"assembly1"],"DE_",grid[i,"pheno"])
    colOfInterest1=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly1"])
    colOfInterest2=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly2"])
    all_results[,outCol]=""
    all_results[!is.na(all_results[,colOfInterest1]) & !is.na(all_results[,colOfInterest2]),outCol]="notSignif"
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest2]<0.05),outCol]="bothSignif"
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest2]>=0.05),outCol]=paste0(grid[i,"assembly1"],"Signif")
    all_results[which(all_results[,colOfInterest1]>=0.05 & all_results[,colOfInterest2]<0.05),outCol]=paste0(grid[i,"assembly2"],"Signif")
    all_results[all_results[,outCol]=="",outCol]=NA
    if(grid[i,"assembly1"]=="T2T"){
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly2"],"Signif"),paste0(grid[i,"assembly1"],"Signif"),"bothSignif","notSignif"))
    }else{
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly1"],"Signif"),paste0(grid[i,"assembly2"],"Signif"),"bothSignif","notSignif"))
        }

    tmp_results=all_results[,c("gene_id",grid[i,"assemblies_DE_pheno"])]
    tmp_results=tmp_results[!is.na(tmp_results[,grid[i,"assemblies_DE_pheno"]]),]
    gene.map = getgo(unlist(lapply(strsplit(unique(tmp_results$gene_id),".",fixed=TRUE),function(x){x[1]})),'hg19','ensGene') 
    
    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly1"],"_",grid[i,"pheno"])
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly1"],"Signif"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
    
    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly2"],"_",grid[i,"pheno"])
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly2"],"Signif"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
}

files=Sys.glob("*module_enrichments.txt")
allGO=list()
for(i in 1:length(files)){
    cat("\r",i,"\t\t\t")
    allGO[[files[i]]]=fread(files[i],data.table=FALSE)
}

signif=lapply(allGO,function(x){x[x[,"BH"]<0.1,]})
signif=do.call("rbind",signif)

##last tab 
for(i in 1:nrow(grid)){
    outCol=paste0(grid[i,"assembly1"],grid[i,"assembly2"],"DE_",grid[i,"pheno"])
    colOfInterest1=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"-",grid[i,"assembly1"])
    colOfInterest1.2=paste0("logFC|",grid[i,"pheno_fixed"],"-",grid[i,"assembly1"])
    colOfInterest2=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"-",grid[i,"assembly2"])
    colOfInterest2.2=paste0("logFC|",grid[i,"pheno_fixed"],"-",grid[i,"assembly2"])
    all_results[,outCol]=""
    all_results[!is.na(all_results[,colOfInterest1]) | !is.na(all_results[,colOfInterest2]),outCol]="notSignif"
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest1.2]>0 & all_results[,colOfInterest2]<0.05 & all_results[,colOfInterest2.2]>0),outCol]="bothSignifUp"
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest1.2]<0 & all_results[,colOfInterest2]<0.05 & all_results[,colOfInterest2.2]<0),outCol]="bothSignifDown"
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest1.2]>0 & all_results[,colOfInterest2]>=0.05),outCol]=paste0(grid[i,"assembly1"],"SignifUp")
    all_results[which(all_results[,colOfInterest1]<0.05 & all_results[,colOfInterest1.2]<0 & all_results[,colOfInterest2]>=0.05),outCol]=paste0(grid[i,"assembly1"],"SignifDown")
    all_results[which(all_results[,colOfInterest1]>=0.05 & all_results[,colOfInterest2]<0.05 & all_results[,colOfInterest2.2]>0),outCol]=paste0(grid[i,"assembly2"],"SignifUp")
    all_results[which(all_results[,colOfInterest1]>=0.05 & all_results[,colOfInterest2]<0.05 & all_results[,colOfInterest2.2]<0),outCol]=paste0(grid[i,"assembly2"],"SignifDown")
    all_results[all_results[,outCol]=="",outCol]=NA
    if(grid[i,"assembly1"]=="T2T"){
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly2"],"SignifUp"),paste0(grid[i,"assembly2"],"SignifDown"),paste0(grid[i,"assembly1"],"SignifUp"),paste0(grid[i,"assembly1"],"SignifDown"),"bothSignifUp","bothSignifDown","notSignif"))
    }else{
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly1"],"SignifUp"),paste0(grid[i,"assembly1"],"SignifDown"),paste0(grid[i,"assembly2"],"SignifUp"),paste0(grid[i,"assembly2"],"SignifDown"),"bothSignifUp","bothSignifDown","notSignif"))
        }

    tmp_results=all_results[,c("gene_id",grid[i,"assemblies_DE_pheno"])]
    tmp_results=tmp_results[!is.na(tmp_results[,grid[i,"assemblies_DE_pheno"]]),]
    gene.map = getgo(unlist(lapply(strsplit(unique(tmp_results$gene_id),".",fixed=TRUE),function(x){x[1]})),'hg19','ensGene') 
    
    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly1"],"_",grid[i,"pheno"],"_Up")
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly1"],"SignifUp"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))

    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly1"],"_",grid[i,"pheno"],"_Down")
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly1"],"SignifDown"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
    
    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly2"],"_",grid[i,"pheno"],"Up")
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly2"],"SignifUp"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))

    name=paste0("DE_AD_",grid[i,"assembly1"],grid[i,"assembly2"],"_",grid[i,"assembly2"],"_",grid[i,"pheno"],"Down")
    tmp_results$DE_status <- ifelse(tmp_results[,grid[i,"assemblies_DE_pheno"]] == paste0(grid[i,"assembly2"],"SignifDown"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
}

grid2=assembliesComb
phenos2=data.frame(pheno=as.character(phenos))
phenos2$pheno_fixed=as.character(phenos2$pheno)
phenos2$pheno_fixed[phenos2$pheno=="ceradsc_defvsctl"]="ceradsc_test.txt"

for(i in 1:nrow(grid2)){
    outCol=paste0(grid2[i,"assembly1"],grid2[i,"assembly2"],"DE_allPhenos")
    colOfInterest1.1=paste0("adj.P.Val|",phenos2[1,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest2.1=paste0("adj.P.Val|",phenos2[2,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest3.1=paste0("adj.P.Val|",phenos2[3,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest4.1=paste0("adj.P.Val|",phenos2[4,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest5.1=paste0("adj.P.Val|",phenos2[5,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest6.1=paste0("adj.P.Val|",phenos2[6,"pheno_fixed"],"-",grid2[i,"assembly1"])
    colOfInterest1.2=paste0("adj.P.Val|",phenos2[1,"pheno_fixed"],"-",grid2[i,"assembly2"])
    colOfInterest2.2=paste0("adj.P.Val|",phenos2[2,"pheno_fixed"],"-",grid2[i,"assembly2"])
    colOfInterest3.2=paste0("adj.P.Val|",phenos2[3,"pheno_fixed"],"-",grid2[i,"assembly2"])
    colOfInterest4.2=paste0("adj.P.Val|",phenos2[4,"pheno_fixed"],"-",grid2[i,"assembly2"])
    colOfInterest5.2=paste0("adj.P.Val|",phenos2[5,"pheno_fixed"],"-",grid2[i,"assembly2"])
    colOfInterest6.2=paste0("adj.P.Val|",phenos2[6,"pheno_fixed"],"-",grid2[i,"assembly2"])
    all_results[,outCol]=""
    all_results[(!is.na(all_results[,colOfInterest1.1]) & !is.na(all_results[,colOfInterest2.1]) & !is.na(all_results[,colOfInterest3.1]) & !is.na(all_results[,colOfInterest4.1]) & !is.na(all_results[,colOfInterest5.1]) & !is.na(all_results[,colOfInterest6.1])) | (!is.na(all_results[,colOfInterest1.2]) & !is.na(all_results[,colOfInterest2.2]) & !is.na(all_results[,colOfInterest3.2]) & !is.na(all_results[,colOfInterest4.2]) & !is.na(all_results[,colOfInterest5.2]) & !is.na(all_results[,colOfInterest6.2])),outCol]="notSignif"
    all_results[which((all_results[,colOfInterest1.1]<0.05 | all_results[,colOfInterest2.1]<0.05 | all_results[,colOfInterest3.1]<0.05 | all_results[,colOfInterest4.1]<0.05 | all_results[,colOfInterest5.1]<0.05 | all_results[,colOfInterest6.1]<0.05) & (all_results[,colOfInterest1.2]<0.05 | all_results[,colOfInterest2.2]<0.05 | all_results[,colOfInterest3.2]<0.05 | all_results[,colOfInterest4.2]<0.05 | all_results[,colOfInterest5.2]<0.05 | all_results[,colOfInterest6.2]<0.05)),outCol]="bothSignif"
    all_results[which((all_results[,colOfInterest1.1]<0.05 | all_results[,colOfInterest2.1]<0.05 | all_results[,colOfInterest3.1]<0.05 | all_results[,colOfInterest4.1]<0.05 | all_results[,colOfInterest5.1]<0.05 | all_results[,colOfInterest6.1]<0.05) & (all_results[,colOfInterest1.2]>=0.05 & all_results[,colOfInterest2.2]>=0.05 & all_results[,colOfInterest3.2]>=0.05 & all_results[,colOfInterest4.2]>=0.05 & all_results[,colOfInterest5.2]>=0.05 & all_results[,colOfInterest6.2]>=0.05)),outCol]=paste0(grid2[i,"assembly1"],"Signif")
    all_results[which((all_results[,colOfInterest1.1]>=0.05 & all_results[,colOfInterest2.1]>=0.05 & all_results[,colOfInterest3.1]>=0.05 & all_results[,colOfInterest4.1]>=0.05 & all_results[,colOfInterest5.1]>=0.05 & all_results[,colOfInterest6.1]>=0.05) & (all_results[,colOfInterest1.2]<0.05 | all_results[,colOfInterest2.2]<0.05 | all_results[,colOfInterest3.2]<0.05 | all_results[,colOfInterest4.2]<0.05 | all_results[,colOfInterest5.2]<0.05 | all_results[,colOfInterest6.2]<0.05)),outCol]=paste0(grid2[i,"assembly2"],"Signif")
    all_results[all_results[,outCol]=="",outCol]=NA
    if(grid2[i,"assembly1"]=="T2T"){
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid2[i,"assembly2"],"Signif"),paste0(grid2[i,"assembly1"],"Signif"),"bothSignif","notSignif"))
    }else{
        all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid2[i,"assembly1"],"Signif"),paste0(grid2[i,"assembly2"],"Signif"),"bothSignif","notSignif"))
        }
        
    tmp_results=all_results[,c("gene_id",outCol)]
    tmp_results=tmp_results[!is.na(tmp_results[,outCol]),]
    gene.map = getgo(unlist(lapply(strsplit(unique(tmp_results$gene_id),".",fixed=TRUE),function(x){x[1]})),'hg19','ensGene') 
    
    name=paste0("DE_AD_",grid2[i,"assembly1"],grid2[i,"assembly2"],"_",grid2[i,"assembly1"],"_allPhenos")
    tmp_results$DE_status <- ifelse(tmp_results[,outCol] == paste0(grid2[i,"assembly1"],"Signif"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
    
    name=paste0("DE_AD_",grid2[i,"assembly1"],grid2[i,"assembly2"],"_",grid2[i,"assembly2"],"_allPhenos")
    tmp_results$DE_status <- ifelse(tmp_results[,outCol] == paste0(grid2[i,"assembly2"],"Signif"), 1, 0)
    df=tmp_results[,c("gene_id", "DE_status")]
    try(annotate_GOterms_DE(mel=df,sample_name=name,gene.map=gene.map,rrvgo=TRUE,revigo=FALSE,GOFigure=FALSE,useFDR=TRUE))
}


files=Sys.glob("*module_enrichments.txt")
allGO=list()
for(i in 1:length(files)){
    cat("\r",i,"\t\t\t")
    allGO[[files[i]]]=fread(files[i],data.table=FALSE)
}

signif=lapply(allGO,function(x){x[x[,"BH"]<0.05,]})
signif=do.call("rbind",signif)
