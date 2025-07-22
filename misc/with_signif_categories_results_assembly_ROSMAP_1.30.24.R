# -----------------------------
# Clear workspace
# -----------------------------
rm(list = ls())  # Remove all existing objects for a clean environment

# -----------------------------
# Load libraries
# -----------------------------
library(data.table)

# -----------------------------
# Load mapping data
# -----------------------------
# mapping: combined mapping between genome assemblies
mapping <- readRDS("rbind_map_between_assemblies_10.04.23.RDS") # find code to generate AD-DEGs-genome-updates/misc/rbind_map_between_assemblies_10.04.23.RDS

# -----------------------------
# Define phenotypes and assemblies
# -----------------------------
# phenos: list of phenotype variables to process
phenos <- c(
  "braaksc_simplified",
  "ceradsc_defvsctl",
  "cogdx_simplified",
  "cts_mmse30"
)
# assemblies: genome assemblies to compare
assemblies <- c("hg18", "hg19", "v30", "v43", "T2T")

# -----------------------------
# Build assembly comparison grid
# -----------------------------
# Create all pairwise combinations of assemblies
assembliesComb <- t(combn(assemblies, 2))
# Separate out those involving T2T for ordering
assembliesComb2 <- assembliesComb[assemblysComb[,2] == "T2T", ]
assembliesComb <- assembliesComb[assemblysComb[,2] != "T2T", ]
# Reverse column order for non-T2T, then append T2T pairs
assembliesComb <- rbind(assembliesComb[, 2:1], assembliesComb2[, 2:1])
# Add combo labels without and with underscore
assembliesComb <- cbind(
  assembliesComb,
  apply(assembliesComb, 1, function(x) paste0(x[1], x[2])),
  apply(assembliesComb, 1, function(x) paste0(x[1], "_", x[2]))
)
colnames(assembliesComb) <- c("assembly1", "assembly2", "combo", "combo_map")

# -----------------------------
# Build grid of assembly-phenotype combinations
# -----------------------------
grid <- expand.grid(
  assemblies = assembliesComb[, "combo"],
  pheno = phenos,
  stringsAsFactors = FALSE
)
# Create column names for DE files and fixed pheno filenames
grid$assemblies_DE_pheno <- paste0(grid$assemblies, "DE_", grid$pheno)
grid$pheno_fixed <- grid$pheno
# Replace one pheno filename to match files on disk
grid$pheno_fixed[grid$pheno_fixed == "ceradsc_defvsctl"] <- "ceradsc_test.txt"
# Merge with assembliesComb to get paths and labels
grid <- merge(grid, assembliesComb, by.x = "assemblies", by.y = "combo")
# Construct path to mapping directory
grid$assembly_map_path <- paste0(
  "path_", grid$assembly1, "_", grid$assembly2, "/"
)

# -----------------------------
# Unique map paths
# -----------------------------
grid2 <- unique(
  grid[, c("assemblies", "combo_map", "assembly_map_path")]
)

# -----------------------------
# Load and merge DE files for each map path
# -----------------------------
for (i in seq_len(nrow(grid2))) {
  path <- grid2[i, "assembly_map_path"]        # Path to DE files
  combo_label <- grid2[i, "combo_map"]         # Label for this assembly pair

  # List all DE files in directory
  DE_files <- list.files(path, pattern = "^DE")
  all <- data.table()  # Initialize empty data.table

  # Loop through each DE file
  for (DE_file in DE_files) {
    tmp <- fread(file.path(path, DE_file))
    if (nrow(all) < 1) {
      all <- tmp  # First file initializes 'all'
    } else {
      all <- merge(all, tmp, by = "V1", all = TRUE)
    }
  }

  # Rename columns to include combo_map label
  new_colnames <- paste0(colnames(all), "_", combo_label)
  colnames(all) <- new_colnames
  colnames(all)[1] <- "gene_id"  # First column is gene ID

  # Assign merged table to variable named by combo label
  assign(combo_label, all)
}

# -----------------------------
# Merge all per-pair results into one data.frame
# -----------------------------
# Start with first pair and iteratively merge all tables
all_results <- Reduce(
  function(x, y) merge(x, get(y), by = "gene_id", all = TRUE),
  as.list(grid2$combo_map)
)
all_results <- as.data.frame(all_results)

# Rename results df and cleanup
assembly_rosmap <- all_results

# -----------------------------
# Define significance status across comparisons
# -----------------------------
# Rebuild assembliesComb with additional labels
assembliesComb <- as.data.frame(
  cbind(
    assembliesComb,
    p_val = "adj.P.Val",
    combo_map2 = paste0("adj.P.Val_", assembliesComb$combo_map)
  ),
  stringsAsFactors = FALSE
)

# Loop to assign significance and common-status flags
for (i in seq_len(nrow(assembliesComb))) {
  combo_label <- assembliesComb$combo_map2[i]
  status_col <- paste0(assembliesComb$combo[i], "DE_status")
  common_col <- paste0(assembliesComb$combo[i], "_common_status")

  # Initialize status columns
  assembly_rosmap[[status_col]] <- ""
  assembly_rosmap[[common_col]] <- TRUE

  # Mark non-significant and significant
  assembly_rosmap[[status_col]][
    is.na(assembly_rosmap[[combo_label]])
  ] <- "notSignif"
  assembly_rosmap[[status_col]][
    assembly_rosmap[[combo_label]] <= 0.05
  ] <- "Signif"

  # Subset mapping to determine common genes
  map_subset <- mapping[mapping$combo == assembliesComb$combo[i], ]
  noncommon <- map_subset[rowSums(is.na(map_subset)) > 0, "common_name"]
  assembly_rosmap[[common_col]][
    assembly_rosmap$gene_id %in% noncommon
  ] <- FALSE
}

# -----------------------------
# Save final results
# -----------------------------
saveRDS(assembly_rosmap,"with_signif_categories_results_assembly_ROSMAP_1.30.24.RDS")
