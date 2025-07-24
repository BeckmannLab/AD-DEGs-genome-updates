library(ggplot2)

matrixhg18 = readRDS("hg18_for_wilcox.RDS") #see AD-DEGs-genome-updates/misc/ref_for_wilcox.R
matrixhg19 = readRDS("hg19_for_wilcox.RDS") #see AD-DEGs-genome-updates/misc/ref_for_wilcox.R
matrixv30 = readRDS("v30_for_wilcox.RDS") #see AD-DEGs-genome-updates/misc/ref_for_wilcox.R
matrixv43 = readRDS("v43_for_wilcox.RDS") #see AD-DEGs-genome-updates/misc/ref_for_wilcox.R
matrixT2T = readRDS("T2T_for_wilcox.RDS") #see AD-DEGs-genome-updates/misc/ref_for_wilcox.R


mats <- list(hg18 = matrixhg18, hg19 = matrixhg19, v30 = matrixv30, v43 = matrixv43, T2T = matrixT2T)

## Functions ------------------------------------------------------------------
align_pair <- function(m1, m2) {
  g <- rownames(m1)[rownames(m1) %in% rownames(m2)]          # genes
  s <- colnames(m1)[colnames(m1) %in% colnames(m2)]          # samples
  list(a = m1[g, s, drop = FALSE],
       b = m2[g, s, drop = FALSE])
}

run_pair <- function(m1, m2, lbl) {
  p   <- align_pair(m1, m2)
  tst <- perform_tests(p$a, p$b)                 # list of htest objects
  pv  <- sapply(tst, `[[`, "p.value")
  vv  <- sapply(tst, `[[`, "statistic")

  out <- data.frame(
    Gene        = names(tst),
    P_Value     = pv,
    V_Value     = vv,
    adj.p.value = p.adjust(pv, "bonferroni"),
    stringsAsFactors = FALSE
  )
  out$Significant <- ifelse(out$adj.p.value <= 0.05, "Yes", "No")
  out$assembly    <- lbl
  out
}

perform_tests <- function(mat1, mat2) {
  stopifnot(identical(dim(mat1), dim(mat2)))  # safety check
  gene_names <- rownames(mat1)
  
  res <- vector("list", length(gene_names))
  for (i in seq_along(gene_names)) {
    res[[i]] <- wilcox.test(
      mat1[i, ], mat2[i, ],
      paired = TRUE, exact = FALSE
    )
  }
  names(res) <- gene_names
  res
}

## Comparison plan ------------------------------------------------------------
pairs <- data.frame(
  a1    = c("hg18","hg18","hg18","hg18",
            "hg19","hg19","hg19",
            "v30","v30",
            "v43"),
  a2    = c("hg19","v30","v43","T2T",
            "v30","v43","T2T",
            "v43","T2T",
            "T2T"),
  label = c(
    "GRCh37/NCBI36",
    "GRCh38.p12/NCBI36",
    "GRCh38.p13/NCBI36",
    "T2T-CHM13v2.0/NCBI36",
    "GRCh38.p12/GRCh37",
    "GRCh38.p13/GRCh37",
    "CHM13v2.0/GRCh37",
    "GRCh38.p13/GRCh38.p12",
    "CHM13v2.0/GRCh38.p12",
    "CHM13v2.0/GRCh38.p13"
  ),
  stringsAsFactors = FALSE
)

## Run tests ------------------------------------------------------------------
results_list <- vector("list", nrow(pairs))

for (i in seq_len(nrow(pairs))) {
	print(i)
  results_list[[i]] <- run_pair(
    mats[[ pairs$a1[i] ]],
    mats[[ pairs$a2[i] ]],
    pairs$label[i]
  )
}

results <- do.call(rbind, results_list)
results <- results[complete.cases(results), ]      # drop rows with any NA

## Plot -----------------------------------------------------------------------
order_vec <- c(
  "GRCh37/NCBI36",
  "GRCh38.p12/NCBI36",
  "GRCh38.p13/NCBI36",
  "T2T-CHM13v2.0/NCBI36",
  "GRCh38.p12/GRCh37",
  "GRCh38.p13/GRCh37",
  "CHM13v2.0/GRCh37",
  "GRCh38.p13/GRCh38.p12",
  "CHM13v2.0/GRCh38.p12",
  "CHM13v2.0/GRCh38.p13"
)

results$assembly <- factor(results$assembly, levels = order_vec)
g <- ggplot(results, aes(assembly, fill = Significant)) +
  geom_bar(position = "dodge") +
  labs(x = "Reference Comparisons", y = "Expressed Gene Count") +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("> 0.05", "<= 0.05"),
                    name   = "Bonferroni") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 50, hjust = 0.8))

ggsave(plot_file, g, useDingbats = FALSE)
