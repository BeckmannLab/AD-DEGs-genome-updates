#Supplementary figure 15
library(tidyr)
library(ggplot2)
library(data.table)

# 1. Define traits and their panel titles
traits <- c(
  "CDR_simplified",
  "CERJ_defvsctl",
  "PlaqueMean_test.txt"
)
panel_titles <- c(
  CDR_simplified       = "CDR",
  CERJ_defvsctl        = "CERAD",
  PlaqueMean_test.txt  = "Plaque mean density"
)

# 2. Read in combo definitions
combo_var <- read.delim("combo_var.txt",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE
)
# structure(list(V1 = 1:10, V2 = c("hg19_hg18", "v30_hg18", "v43_hg18", 
# "T2T_hg18", "v30_hg19", "v43_hg19", "T2T_hg19", "v30_v43", "T2T_v30", 
# "T2T_v43")), class = "data.frame", row.names = c(NA, -10L))

combo_var[,1] <- c(
  "GRCh37/NCBI36","GRCh38.p12/NCBI36","GRCh38.p13/NCBI36","CHM13v2.0/NCBI36",
  "GRCh38.p12/GRCh37","GRCh38.p13/GRCh37","CHM13v2.0/GRCh37",
  "GRCh38.p13/GRCh38.p12","CHM13v2.0/GRCh38.p12","CHM13v2.0/GRCh38.p13"
)

# 3. Read in the var lookup table and standardize names
var <- read.delim("var.txt",
  header = FALSE, sep = "\t", stringsAsFactors = FALSE
)

# structure(list(V1 = c("hg18hg19DE_ceradsc_defvsctl", "hg18v30DE_ceradsc_defvsctl", 
# "hg18v43DE_ceradsc_defvsctl", "hg19v30DE_ceradsc_defvsctl", "hg19v43DE_ceradsc_defvsctl", 
# "v30v43DE_ceradsc_defvsctl", "T2Thg18DE_ceradsc_defvsctl", "T2Thg19DE_ceradsc_defvsctl", 
# "T2Tv30DE_ceradsc_defvsctl", "T2Tv43DE_ceradsc_defvsctl"), V2 = c("hg18Signif", 
# "hg18Signif", "hg18Signif", "hg19Signif", "hg19Signif", "v30Signif", 
# "hg18Signif", "hg19Signif", "v30Signif", "v43Signif"), V3 = c("hg19Signif", 
# "v30Signif", "v43Signif", "v30Signif", "v43Signif", "v43Signif", 
# "T2TSignif", "T2TSignif", "T2TSignif", "T2TSignif")), class = "data.frame", row.names = c(NA, 
# -10L))

var$V1 <- gsub("hg18hg19DE",  "hg19hg18DE", var$V1)
var$V1 <- gsub("hg18v30DE",   "v30hg18DE",  var$V1)
var$V1 <- gsub("hg18v43DE",   "v43hg18DE",  var$V1)
var$V1 <- gsub("hg19v30DE",   "v30hg19DE",  var$V1)
var$V1 <- gsub("hg19v43DE",   "v43hg19DE",  var$V1)
var$V1 <- gsub("v30v43DE",    "v43v30DE",   var$V1)

############################
####DE AD results first
################################

# Define file paths for differential expression results for each genome assembly
path_hg18= #path to DE results
path_hg19= #path to DE results
path_v30= #path to DE results
path_v43= #path to DE results
path_T2T= #path to DE results

# Specify phenotypes or tests to process
phenos=c("CDR_simplified", "CERJ_defvsctl","PlaqueMean_test.txt")

# List of genome assemblies
assemblies=c("hg18","hg19","v30","v43","T2T")

# Generate all unique pairs of assemblies
assembliesComb = t(combn(assemblies,2))
assembliesComb2=assembliesComb[assembliesComb[,2]=="T2T",]
assembliesComb=assembliesComb[assembliesComb[,2]!="T2T",]
assembliesComb=rbind(assembliesComb,assembliesComb2[,2:1])
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[1],x[2])}))

# Create a grid of assembly pairs and phenotypes
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")

# Prepare paths for each assembly
grid3=as.data.frame(assemblies)
grid3$assemblies_path=paste("path",grid3$assemblies, sep = "_")

# Loop through each assembly to read and merge DE files
for(i in 1:nrow(grid3)){
  path=get(grid3[i,"assemblies_path"])
  DE_files=list.files(path,pattern="^DE")
  all=data.table()
  # Read and merge each DE file
  for(DE_file in DE_files){
      tmp=fread(paste0(path,DE_file))
      phenotype=unlist(lapply(strsplit(DE_file,"_"),function(x){paste0(x[11],"_",x[12])}))
      colnames(tmp)[2:ncol(tmp)]=paste0(colnames(tmp)[2:ncol(tmp)],"|",phenotype)
      if(nrow(all)<1){
          all=tmp
      }else{
          all=merge(all,tmp,by="V1",all=TRUE)
      }}
  colnames(all)=paste0(colnames(all),"_",grid3[i,"assemblies"])
  colnames(all)[1]="gene_id"
  value_assigned=as.character(grid3[i, "assemblies"])
  assign(value_assigned,all)
}

# Merge results from all assemblies into one data frame
all_results_AD=merge(hg18,hg19,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,v30,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,v43,by="gene_id",all=TRUE)
all_results_AD=merge(all_results_AD,T2T,by="gene_id",all=TRUE)
all_results_AD=as.data.frame(all_results_AD)

# Reinitialize combinations for significance categorization
assembliesComb = t(combn(assemblies,2))
assembliesComb=cbind(assembliesComb,apply(assembliesComb,1,function(x){paste0(x[2],x[1])}))

# Recreate grid for significance loop
grid=expand.grid(assembliesComb[,3],phenos)
grid$Var3=paste0(grid$Var1,"DE_",grid$Var2)
grid$Var4=as.character(grid$Var2)
colnames(grid)=c("assemblies","pheno","assemblies_DE_pheno","pheno_fixed")
colnames(assembliesComb)=c("assembly1","assembly2","assemblies")
grid=merge(grid,assembliesComb,by="assemblies")

# Copy merged results for analysis
all_results = all_results_AD

# Loop to categorize significance across assembly pairs for each phenotype
for(i in 1:nrow(grid)){
    # Define output column for this comparison
    outCol=paste0(grid[i,"assembly2"],grid[i,"assembly1"],"DE_",grid[i,"pheno"])
    # Identify adjusted p-value columns for each assembly
    colOfInterest1=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly1"])
    colOfInterest2=paste0("adj.P.Val|",grid[i,"pheno_fixed"],"_",grid[i,"assembly2"])
    # Initialize output column
    all_results[,outCol]=""

    # Default label when both present but not significant
    all_results[!is.na(all_results[,colOfInterest1]) & !is.na(all_results[,colOfInterest2]),outCol]="notSignif"

    # Cases where one is NA and the other significant
    all_results[which(is.na(all_results[, colOfInterest1]) & all_results[, colOfInterest2] <= 0.05), outCol] = paste0(grid[i,"assembly2"],"Signif")
    all_results[which(is.na(all_results[,colOfInterest2]) & all_results[,colOfInterest1]<=0.05),outCol]=paste0(grid[i,"assembly1"],"Signif")

    # Both significant
    all_results[which(all_results[,colOfInterest1]<=0.05 & all_results[,colOfInterest2]<=0.05),outCol]="bothSignif"

    # One significant, one not
    all_results[which(all_results[,colOfInterest1]<=0.05 & all_results[,colOfInterest2]>0.05),outCol]=paste0(grid[i,"assembly1"],"Signif")
    all_results[which(all_results[,colOfInterest1]>0.05 & all_results[,colOfInterest2]<=0.05),outCol]=paste0(grid[i,"assembly2"],"Signif")

    # Cases where one is NA but not significant
    all_results[which(is.na(all_results[,colOfInterest1]) & all_results[,colOfInterest2]>0.05),outCol]="notSignif"
    all_results[which(is.na(all_results[,colOfInterest2]) & all_results[,colOfInterest1]>0.05),outCol]="notSignif"

    # Replace empty strings with NA
    all_results[all_results[,outCol]=="",outCol]=NA

    # Set factor levels for categorization
    all_results[,outCol]=factor(all_results[,outCol],levels=c(paste0(grid[i,"assembly1"],"Signif"),paste0(grid[i,"assembly2"],"Signif"),"bothSignif","notSignif"))
}

# Convert to data frame and save results
all_results=as.data.frame(all_results)

# 5. Build a summary for each trait
overview_list <- list()
for (tr in traits) {
  # Update var$V1 so it points to the correct columns for this trait
  var_tr <- var
  var_tr$V1 <- gsub("ceradsc_defvsctl", tr, var_tr$V1)
  
  # Subset columns that match this trait
  subset_tr <- all_results[, grep(tr, colnames(all_results)), drop = FALSE]
  
  # Initialize a data.frame to hold counts
  df <- data.frame(
    assembly = combo_var[,1],
    older    = integer(nrow(combo_var)),
    newer    = integer(nrow(combo_var)),
    both     = integer(nrow(combo_var)),
    trait    = tr,
    stringsAsFactors = FALSE
  )
  
  # Fill in counts for each assembly
  for (i in seq_len(nrow(var_tr))) {
    vec <- subset_tr[[ var_tr$V1[i] ]]
    df$older[i] <- sum(vec == var_tr$V2[i],    na.rm = TRUE)
    df$newer[i] <- sum(vec == var_tr$V3[i],    na.rm = TRUE)
    df$both[i]  <- sum(vec == "bothSignif",    na.rm = TRUE)
  }
  
  overview_list[[tr]] <- df
}

# 6. Combine and pivot to long format
overview_all <- do.call(rbind, overview_list)
overview_t   <- pivot_longer(
  overview_all,
  cols      = c("older", "newer", "both"),
  names_to  = "comparison",
  values_to = "value"
)

# 7. Factor levels for consistent ordering
overview_t$assembly <- factor(
  overview_t$assembly,
  levels = combo_var[,1]
)
overview_t$trait <- factor(
  overview_t$trait,
  levels = traits
)

# ==== FIX: rename & factor the 'comparison' column instead of 'name' ====
overview_t$comparison <- gsub("older", "Older reference", overview_t$comparison)
overview_t$comparison <- gsub("newer", "Newer reference", overview_t$comparison)
overview_t$comparison <- gsub("both",  "Both references", overview_t$comparison)

overview_t$comparison <- factor(
  overview_t$comparison,
  levels = c("Older reference", "Newer reference", "Both references")
)

# 8. Create the faceted barplot with 3 columns
g <- ggplot(overview_t, aes(x = assembly, y = value, fill = comparison)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    breaks = c("Older reference", "Newer reference", "Both references"),
    values = c(
      "Older reference"  = "darkblue",
      "Newer reference"  = "lightblue",
      "Both references"  = "darkgreen"
    )
  ) +
  labs(
    x    = "Reference Comparison",
    y    = "Number of DEGs",
    fill = "AD DEG in"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 19),
    axis.text         = element_text(size = 22),
    axis.title        = element_text(size = 24),
    legend.title      = element_text(size = 20),
    legend.text       = element_text(size = 18),
    strip.text        = element_text(size = 20),
    legend.background = element_rect(fill = "white", colour = "black", size = 0.5),
    legend.key        = element_rect(fill = NA, colour = NA)
  ) +
  facet_wrap(
    ~ trait,
    ncol     = 3,
    labeller = as_labeller(panel_titles),
    scales = "free_y"
  )

# 9. Save the plot
ggsave(
  filename    = "AD_DEGs_msbb_direction_all_traits.pdf",
  plot        = g,
  useDingbats = FALSE,
  width       = 16,
  height      = 8
)
