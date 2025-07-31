# Load required libraries
library(ggplot2)
options(stringsAsFactors = FALSE)
library("R.matlab")
library(goseq)
library(topGO)
library(org.Hs.eg.db)
library(Rgraphviz)
library(data.table)

# -------------------------------------
# Load data and define parameters
# -------------------------------------

interaction_rosmap <- readRDS("results_interaction_ROSMAP_1.8.24") #see AD-DEGs-genome-updates/misc/results_interaction_ROSMAP_1.8.24.R
mapping <- readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.R

phenos <- c("braaksc_simplified", "ceradsc_defvsctl","cogdx_simplified", "cts_mmse30")
assemblies <- c("hg18", "hg19", "v30", "v43", "T2T")

main_path <- "interaction_allPathways/for_each_pheno/"
setwd(main_path)

# Generate all pairwise combinations of assemblies
assembliesComb <- t(combn(assemblies, 2))
assembliesComb2 <- assembliesComb[assembliesComb[,2] == "T2T",]
assembliesComb <- rbind(assembliesComb[,2:1], assembliesComb2[,2:1])
assembliesComb <- cbind(
  assembliesComb,
  apply(assembliesComb, 1, function(x) paste0(x[1], x[2])),
  apply(assembliesComb, 1, function(x) paste0(x[1], "_", x[2]))
)
colnames(assembliesComb) <- c("assembly1", "assembly2", "assemblies", "assembly_map")

# Build a grid of all (assembly pair x phenotype) combinations
grid <- expand.grid(assembliesComb[, "assemblies"], phenos)
grid$Var3 <- paste0(grid$Var1, "DE_", grid$Var2)
grid$Var4 <- as.character(grid$Var2)
grid$Var4[grid$Var4 == "ceradsc_defvsctl"] <- "ceradsc_test.txt"
colnames(grid) <- c("assemblies", "pheno", "assemblies_DE_pheno", "pheno_fixed")
grid <- merge(grid, assembliesComb, by = "assemblies")

# -------------------------------------
# Run GO enrichment for each phenotype
# -------------------------------------

for (i in 1:nrow(grid)) {
  outCol <- paste0(grid[i,"assembly1"], grid[i,"assembly2"], "DE_", grid[i,"pheno"])
  colOfInterest1 <- paste0("adj.P.Val|", grid[i,"pheno"], "_", grid[i,"assembly1"], "_", grid[i,"assembly2"])
  colOfInterest2 <- paste0("logFC|", grid[i,"pheno"], "_", grid[i,"assembly1"], "_", grid[i,"assembly2"])

  interaction_rosmap[, outCol] <- ""
  interaction_rosmap[!is.na(interaction_rosmap[, colOfInterest1]), outCol] <- "notSignif"
  interaction_rosmap[which(interaction_rosmap[, colOfInterest1] <= 0.05 & interaction_rosmap[, colOfInterest2] > 0), outCol] <- paste0(grid[i,"assembly1"], "Signif")
  interaction_rosmap[which(interaction_rosmap[, colOfInterest1] <= 0.05 & interaction_rosmap[, colOfInterest2] <= 0), outCol] <- paste0(grid[i,"assembly2"], "Signif")
  interaction_rosmap[interaction_rosmap[, outCol] == "", outCol] <- NA

  interaction_rosmap[, outCol] <- factor(interaction_rosmap[, outCol],
    levels = c(paste0(grid[i,"assembly1"],"Signif"), paste0(grid[i,"assembly2"],"Signif"), "notSignif"))

  mapping_subset <- mapping[mapping$combo == grid$assembly_map[i],]
  mapping_subset2 <- mapping_subset[rowSums(is.na(mapping_subset)) > 0,]
  interaction_rosmap[interaction_rosmap$gene_id %in% mapping_subset2$common_name, outCol] <- NA

  tmp_results <- interaction_rosmap[, c("gene_id", grid[i,"assemblies_DE_pheno"])]
  tmp_results <- tmp_results[!is.na(tmp_results[,2]), ]
  gene.map <- getgo(unlist(lapply(strsplit(unique(tmp_results$gene_id), ".", fixed=TRUE), function(x) x[1])), 'hg19', 'ensGene')

  name1 <- paste0("interaction_", grid[i,"assembly1"], grid[i,"assembly2"], "_", grid[i,"assembly1"], "_", grid[i,"pheno"])
  tmp_results$DE_status <- ifelse(tmp_results[,2] == paste0(grid[i,"assembly1"], "Signif"), 1, 0)
  df <- tmp_results[, c("gene_id", "DE_status")]
  try(annotate_GOterms_DE(mel = df, sample_name = name1, gene.map = gene.map, rrvgo = TRUE, revigo = FALSE, GOFigure = FALSE, useFDR = TRUE))

  name2 <- paste0("interaction_", grid[i,"assembly1"], grid[i,"assembly2"], "_", grid[i,"assembly2"], "_", grid[i,"pheno"])
  tmp_results$DE_status <- ifelse(tmp_results[,2] == paste0(grid[i,"assembly2"], "Signif"), 1, 0)
  df <- tmp_results[, c("gene_id", "DE_status")]
  try(annotate_GOterms_DE(mel = df, sample_name = name2, gene.map = gene.map, rrvgo = TRUE, revigo = FALSE, GOFigure = FALSE, useFDR = TRUE))
}

# -------------------------------------
# Aggregate and save per-phenotype results
# -------------------------------------

files <- Sys.glob("*module_enrichments.txt")
allGO <- list()
for (i in 1:length(files)) {
  cat("\r", i, "\t\t\t")
  allGO[[files[i]]] <- fread(files[i], data.table = FALSE)
}

signif <- lapply(allGO, function(x) x[x[,"BH"] < 20, ])
signif <- do.call("rbind", signif)
signif <- signif[order(signif$BH), ]

signif$assembly_pairs <- unlist(lapply(strsplit(rownames(signif), "_"), function(x) x[2]))
signif$phenotype <- unlist(lapply(strsplit(rownames(signif), "_"), function(x) paste(x[4], x[5], sep = "_")))

saveRDS(signif, "signif_rosmap.RDS")


###############
#across traits 
###############


# Load interaction data and mapping
interaction_rosmap <- readRDS("results_interaction_ROSMAP_1.8.24") #see AD-DEGs-genome-updates/misc/results_interaction_ROSMAP_1.8.24.R
mapping <- readRDS("rbind_map_between_assemblies_10.04.23.RDS") #see AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.R

# Define output directory and set working directory
main_path <- "across_all_phenos/"
system(paste("mkdir", main_path))
setwd(main_path)

# Reuse assemblies and phenotypes
grid2 <- assembliesComb
phenos2 <- data.frame(pheno = as.character(phenos))
phenos2$pheno_fixed <- phenos2$pheno
phenos2$pheno_fixed[phenos2$pheno == "ceradsc_defvsctl"] <- "ceradsc_test.txt"

# Loop over each assembly pair across all phenotypes
for (i in 1:nrow(grid2)) {
  outCol <- paste0(grid2[i, "assembly1"], grid2[i, "assembly2"], "DE_allPhenos")

  # Construct column names for each phenotype
  colOfInterest1.1 <- paste0("adj.P.Val|", phenos2[1,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest2.1 <- paste0("adj.P.Val|", phenos2[2,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest3.1 <- paste0("adj.P.Val|", phenos2[3,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest4.1 <- paste0("adj.P.Val|", phenos2[4,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest5.1 <- paste0("adj.P.Val|", phenos2[5,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest6.1 <- paste0("adj.P.Val|", phenos2[6,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])

  colOfInterest1.2 <- paste0("logFC|", phenos2[1,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest2.2 <- paste0("logFC|", phenos2[2,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest3.2 <- paste0("logFC|", phenos2[3,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest4.2 <- paste0("logFC|", phenos2[4,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest5.2 <- paste0("logFC|", phenos2[5,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])
  colOfInterest6.2 <- paste0("logFC|", phenos2[6,"pheno"], "_", grid2[i,"assembly1"], "_", grid2[i,"assembly2"])

  interaction_rosmap[, outCol] <- ""

  # Assign "notSignif" if any adjusted p-values are present
  interaction_rosmap[(!is.na(interaction_rosmap[, colOfInterest1.1]) |
                      !is.na(interaction_rosmap[, colOfInterest2.1]) |
                      !is.na(interaction_rosmap[, colOfInterest3.1]) |
                      !is.na(interaction_rosmap[, colOfInterest4.1]) |
                      !is.na(interaction_rosmap[, colOfInterest5.1]) |
                      !is.na(interaction_rosmap[, colOfInterest6.1])),
                     outCol] <- "notSignif"

  # Mark as assembly1Signif if any phenotype shows positive logFC and p < 0.05
  interaction_rosmap[which((interaction_rosmap[, colOfInterest1.1] < 0.05 & interaction_rosmap[, colOfInterest1.2] > 0) |
                           (interaction_rosmap[, colOfInterest2.1] < 0.05 & interaction_rosmap[, colOfInterest2.2] > 0) |
                           (interaction_rosmap[, colOfInterest3.1] < 0.05 & interaction_rosmap[, colOfInterest3.2] > 0) |
                           (interaction_rosmap[, colOfInterest4.1] < 0.05 & interaction_rosmap[, colOfInterest4.2] > 0) |
                           (interaction_rosmap[, colOfInterest5.1] < 0.05 & interaction_rosmap[, colOfInterest5.2] > 0) |
                           (interaction_rosmap[, colOfInterest6.1] < 0.05 & interaction_rosmap[, colOfInterest6.2] > 0)),
                     outCol] <- paste0(grid2[i, "assembly1"], "Signif")

  # Mark as assembly2Signif if logFC <= 0
  interaction_rosmap[which((interaction_rosmap[, colOfInterest1.1] < 0.05 & interaction_rosmap[, colOfInterest1.2] <= 0) |
                           (interaction_rosmap[, colOfInterest2.1] < 0.05 & interaction_rosmap[, colOfInterest2.2] <= 0) |
                           (interaction_rosmap[, colOfInterest3.1] < 0.05 & interaction_rosmap[, colOfInterest3.2] <= 0) |
                           (interaction_rosmap[, colOfInterest4.1] < 0.05 & interaction_rosmap[, colOfInterest4.2] <= 0) |
                           (interaction_rosmap[, colOfInterest5.1] < 0.05 & interaction_rosmap[, colOfInterest5.2] <= 0) |
                           (interaction_rosmap[, colOfInterest6.1] < 0.05 & interaction_rosmap[, colOfInterest6.2] <= 0)),
                     outCol] <- paste0(grid2[i, "assembly2"], "Signif")

  interaction_rosmap[interaction_rosmap[, outCol] == "", outCol] <- NA

  mapping_subset <- mapping[mapping$combo == grid2$assembly_map[i], ]
  mapping_subset <- mapping_subset[which(is.na(mapping_subset), arr.ind = TRUE), ]
  interaction_rosmap[interaction_rosmap$gene_id %in% mapping_subset$common_name, outCol] <- NA

  interaction_rosmap[, outCol] <- factor(interaction_rosmap[, outCol],
    levels = c(paste0(grid2[i,"assembly1"],"Signif"), paste0(grid2[i,"assembly2"],"Signif"), "notSignif"))

  tmp_results <- interaction_rosmap[, c("gene_id", outCol)]
  tmp_results <- tmp_results[!is.na(tmp_results[, outCol]), ]
  gene.map <- getgo(unlist(lapply(strsplit(unique(tmp_results$gene_id), ".", fixed = TRUE), function(x) x[1])), 'hg19', 'ensGene')

  name1 <- paste0("Interaction_", grid2[i,"assembly1"], grid2[i,"assembly2"], "_", grid2[i,"assembly1"], "_allPhenos")
  tmp_results$DE_status <- ifelse(tmp_results[, outCol] == paste0(grid2[i,"assembly1"], "Signif"), 1, 0)
  df <- tmp_results[, c("gene_id", "DE_status")]
  try(annotate_GOterms_DE(mel = df, sample_name = name1, gene.map = gene.map, rrvgo = TRUE, revigo = FALSE, GOFigure = FALSE, useFDR = TRUE))

  name2 <- paste0("Interaction_", grid2[i,"assembly1"], grid2[i,"assembly2"], "_", grid2[i,"assembly2"], "_allPhenos")
  tmp_results$DE_status <- ifelse(tmp_results[, outCol] == paste0(grid2[i,"assembly2"], "Signif"), 1, 0)
  df <- tmp_results[, c("gene_id", "DE_status")]
  try(annotate_GOterms_DE(mel = df, sample_name = name2, gene.map = gene.map, rrvgo = TRUE, revigo = FALSE, GOFigure = FALSE, useFDR = TRUE))
}

# -------------------------------------
# Aggregate and save across-all-phenotypes results
# -------------------------------------

files <- Sys.glob("*module_enrichments.txt")
allGO <- list()
for (i in 1:length(files)) {
  cat("\r", i, "\t\t\t")
  allGO[[files[i]]] <- fread(files[i], data.table = FALSE)
}

signif <- lapply(allGO, function(x) x[x[,"BH"] < 0.05, ])
signif <- do.call("rbind", signif)
signif <- signif[order(signif$BH), ]

signif$assembly_pairs <- unlist(lapply(strsplit(rownames(signif), "_"), function(x) x[3]))
signif$phenotype <- unlist(lapply(strsplit(rownames(signif), "_"), function(x) paste(x[5], sep = "_")))

saveRDS(signif, "/across_all_phenos/signif_rosmap.RDS")
