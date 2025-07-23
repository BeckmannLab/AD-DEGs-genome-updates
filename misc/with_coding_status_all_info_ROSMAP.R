rm(list=ls())
# Load ggplot2 for plotting, data.table for fast file reading, tidyr for data reshaping,
# GGally for ggpairs, and stringr for string manipulation
library(ggplot2)
library(data.table)
library(tidyr)
library(GGally)
library(stringr)

## Differential expression (DE) outputs - direct from toptable
# Read DE results for various assembly comparisons into separate data.frames
hg18_hg19_DE = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_v30_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_v30_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v30_v43_DE   = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_v43_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_v43_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg18_T2T_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
hg19_T2T_DE  = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v30_T2T_DE   = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")
v43_T2T_DE   = read.delim("DE_assembly_corrected_for_id_assembly_test.txt")

## GTF annotation files for each assembly
hg18_gtf = fread("Homo_sapiens.NCBI36.53.gtf")  # NCBI36/hg18 annotation
hg19_gtf = fread("Homo_sapiens.GRCh37.70.gtf")  # GRCh37/hg19 annotation
v30_gtf  = fread("gencode.v30.primary_assembly.annotation.gtf")
v43_gtf  = fread("gencode.v43.primary_assembly.annotation.gtf", data.table=FALSE)
t2t_gtf  = fread("Homo_sapiens_GCA_009914755.4_2022_07_genes.gtf", data.table=FALSE)

## Load mapping between assemblies for common gene names
mapping = readRDS("rbind_map_between_assemblies_10.04.23.RDS")  # map file created via misc script

### Extract coding status per gene for each assembly ###

# hg18
gt = hg18_gtf\ ncolnames(gt) = c("seqname","source","region","start","end","score","strand","frame","V9")
# split V9 on ';' to new columns
attrs = data.frame(do.call(rbind, strsplit(as.character(gt$V9), ";", fixed=TRUE)))
gt2 = cbind(gt, attrs)
# select gene_id and type column, then split gene_id field
gt_sub = gt2[, c("X1","source")]
gt_id = data.frame(do.call(rbind, strsplit(as.character(gt_sub$X1), " ", fixed=TRUE)))
gt_clean = cbind(gene_id_raw=gt_id$V2, status=gt_sub$source)
# remove quotes
gt_clean = transform(gene_id_raw = gsub('"','', gt_clean$gene_id_raw))
hg18_final = data.frame(gene_id = gsub('"', '', gt_clean$gene_id_raw), status = gt_clean$status, stringsAsFactors=FALSE)
# collapse non-protein_coding entries
hg18_final$status[hg18_final$status != "protein_coding"] <- "non_coding"
hg18_final = unique(hg18_final)

# hg19
gt = hg19_gtf
colnames(gt) = c("seqname","source","region","start","end","score","strand","frame","V9")
attrs = data.frame(do.call(rbind, strsplit(as.character(gt$V9), ";", fixed=TRUE)))
gt2 = cbind(gt, attrs)
gt2 = gt2[gt2$region=="exon",]
# extract gene_id and gene_type fields
id_field = str_replace_all(str_extract_all(gt2$X1, '"([^"]+)"'), '"','')
type_field = str_replace_all(str_extract_all(gt2$X5, '"([^"]+)"'), '"','')
hg19_df = data.frame(gene_id=id_field, status=type_field, stringsAsFactors=FALSE)
hg19_df$status[hg19_df$status != "protein_coding"] <- "non_coding"
hg19_final = unique(hg19_df)

# v30
gt = v30_gtf
colnames(gt) = c("seqname","source","region","start","end","score","strand","frame","V9")
gt2 = gt[gt$region=="exon",]
gene_type = sub(".*gene_type\\s(.*?);.*", "\\1", gt2$V9)
gene_id   = sub(".*gene_id\\s(.*?);.*",   "\\1", gt2$V9)
df30 = data.frame(gene_id = gsub('"','',gene_id), status = gsub('"','',gene_type), stringsAsFactors=FALSE)
df30$status[df30$status != "protein_coding"] <- "non_coding"
v30_tmp = unique(df30)
# strip version suffix (e.g. ENSG000001234.5 -> ENSG000001234)
v30_tmp$gene_id = sapply(strsplit(v30_tmp$gene_id, ".", fixed=TRUE), `[`, 1)
v30_final = v30_tmp

# v43
gt = v43_gtf
colnames(gt) = c("seqname","source","region","start","end","score","strand","frame","V9")
gt2 = gt[gt$region=="exon",]
gene_type = sub(".*gene_type\\s(.*?);.*", "\\1", gt2$V9)
gene_id   = sub(".*gene_id\\s(.*?);.*",   "\\1", gt2$V9)
df43 = data.frame(gene_id = gsub('"','',gene_id), status = gsub('"','',gene_type), stringsAsFactors=FALSE)
df43$status[df43$status != "protein_coding"] <- "non_coding"
v43_tmp = unique(df43)
# strip version
v43_tmp$gene_id = sapply(strsplit(v43_tmp$gene_id, ".", fixed=TRUE), `[`, 1)
v43_final = v43_tmp

# T2T: extract gene_biotype instead of gene_type
gt = t2t_gtf
colnames(gt) = c("seqname","source","region","start","end","score","strand","frame","V9")
gt2 = gt[gt$region=="exon",]
gene_biotype = sub(".*gene_biotype\\s(.*?);.*","\\1", gt2$V9)
gene_id      = sub(".*gene_id\\s(.*?);.*","\\1", gt2$V9)
dfT2T = data.frame(gene_id=gsub('"','',gene_id), status=gsub('"','',gene_biotype), stringsAsFactors=FALSE)
dfT2T$status[dfT2T$status != "protein_coding"] <- "non_coding"
t2t_final = unique(dfT2T)

## Prepare comparison grids and merge DE with coding status
# Non-T2T
pairs      = c("hg18_hg19","hg18_v30","hg19_v30","v30_v43","hg18_v43","hg19_v43")
de_out     = c("hg18_hg19_DE","hg18_v30_DE","hg19_v30_DE","v30_v43_DE","hg18_v43_DE","hg19_v43_DE")
combos     = c("GRCh37-NCBI36","GRCh38.12-NCBI36","GRCh38.12-GRCh37","GRCh38.13-GRCh38.12","GRCh38.13-NCBI36","GRCh38.13-GRCh37")
gtf1      = list(hg19_final, v30_final, v30_final, v43_final, v43_final, v43_final)
gtf2      = list(hg18_final, hg18_final, hg19_final, v30_final, hg18_final, hg19_final)

grid = data.frame(pairs, de_out, combos, stringsAsFactors=FALSE)
final_non_t2t = vector("list", nrow(grid))
for(i in seq_len(nrow(grid))) {
  message("Processing ", grid$pairs[i])
  de = get(grid$de_out[i])  # DE table
  # map common vs non-common genes
  m    = mapping[grep(grid$pairs[i], mapping$combo),]
  common   = na.omit(m)
  not_common = m[!complete.cases(m),]
  # annotate DE
  de$sig_status = ifelse(de$adj.P.Val < 0.05, "Signif","notSignif")
  de$direction  = ifelse(de$logFC>0, "New", ifelse(de$logFC<0, "Old", NA))
  de$in_common  = !de$X %in% not_common$common_name
  de$assembly   = grid$combos[i]
  colnames(de)[1] = "gene_id"
  de_sub = de[,c("gene_id","sig_status","direction","in_common","assembly")]
  # merge with coding status tables
  ref_all   = merge(gtf1[[i]], gtf2[[i]], by="gene_id", all=TRUE)
  ref_de    = merge(ref_all, de_sub, by="gene_id", all=TRUE)
  final_non_t2t[[i]] = ref_de
}
combined_df_non_T2T = do.call(rbind, final_non_t2t)

# T2T comparisons
grid_t2t = data.frame(
  pairs = c("hg18_T2T","hg19_T2T","v30_T2T","v43_T2T"),
  de_out = c("hg18_T2T_DE","hg19_T2T_DE","v30_T2T_DE","v43_T2T_DE"),
  combos = c("CHM13v2.0-NCBI36","CHM13v2.0-GRCh37","CHM13v2.0-GRCh38.12","CHM13v2.0-GRCh38.13"),
  stringsAsFactors=FALSE
)
final_t2t = vector("list", nrow(grid_t2t))
for(i in seq_len(nrow(grid_t2t))) {
  message("Processing ", grid_t2t$pairs[i])
  de = get(grid_t2t$de_out[i])
  m    = mapping[grep(grid_t2t$pairs[i], mapping$combo),]
  common   = na.omit(m)
  not_common = m[!complete.cases(m),]
  de$sig_status = ifelse(de$adj.P.Val<0.05,"Signif","notSignif")
  de$direction  = ifelse(de$logFC>0,"New", ifelse(de$logFC<0,"Old",NA))
  de$in_common  = !de$X %in% not_common$common_name
  de$assembly   = grid_t2t$combos[i]
  colnames(de)[1] = "gene_id"
  de_sub = de[,c("gene_id","sig_status","direction","in_common","assembly")]
  # merge T2T status via mapping
  t2t_map = mapping[mapping$combo==grid_t2t$pairs[i],]
  t2t_status = merge(t2t_final, t2t_map, by.x="gene_id", by.y="assembly_2_orginal_gene_name", all=TRUE)
  t2t_status = unique(data.frame(gene_id = t2t_status$common_name, status = t2t_status$status))
  # combine reference and DE
  ref_base = get(ifelse(grepl("hg18",grid_t2t$pairs[i]),"hg18_final", 
                ifelse(grepl("hg19",grid_t2t$pairs[i]),"hg19_final", 
                ifelse(grepl("v30",grid_t2t$pairs[i]),"v30_final","v43_final"))))
  ref_comb  = merge(ref_base, t2t_status, by="gene_id", all=TRUE)
  ref_de    = merge(ref_comb, de_sub, by="gene_id", all=TRUE)
  final_t2t[[i]] = ref_de
}
combined_df_T2T = do.call(rbind, final_t2t)

# Merge all
all_info_ROSMAP = rbind(combined_df_non_T2T, combined_df_T2T)
# filter complete cases for plotting
df_clean = all_info_ROSMAP[complete.cases(all_info_ROSMAP[,c("in_common","assembly","direction","sig_status","status")]),]
# build combined direction-status-in_common field
all_info_ROSMAP = df_clean
after_filter = with(all_info_ROSMAP,
  data.frame(
    direction_signif = ifelse(sig_status=="Signif", direction, "notSignif"),
    in_common       = in_common,
    assembly        = assembly,
    status          = status,
    stringsAsFactors=FALSE
  )
)
all_info_ROSMAP = transform(all_info_ROSMAP,
  for_plot = paste0(direction_signif, "_", in_common)
)
# set factor levels for plotting order
all_info_ROSMAP$assembly = factor(all_info_ROSMAP$assembly,
  levels = c(
    "GRCh37-NCBI36","GRCh38.12-NCBI36","GRCh38.13-NCBI36","CHM13v2.0-NCBI36",
    "GRCh38.12-GRCh37","GRCh38.13-GRCh37","CHM13v2.0-GRCh37",
    "GRCh38.13-GRCh38.12","CHM13v2.0-GRCh38.12","CHM13v2.0-GRCh38.13"
  )
)

## Save final annotated data for downstream plotting or analysis
saveRDS(all_info_ROSMAP, "with_coding_status_all_info_ROSMAP.RDS")
