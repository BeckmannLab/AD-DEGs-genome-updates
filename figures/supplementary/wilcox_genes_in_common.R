# ============================ #
#       LOAD LIBRARIES        #
# ============================ #
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(edgeR)
  library(variancePartition)
  library(ggplot2)
  library(assertthat)
  library(tidyverse)
  library(dplyr)
})

# ============================ #
#         FUNCTIONS            #
# ============================ #

add_decimal_numbers <- function(target_values, value_list) {
  target_indices <- match(value_list, target_values)
  ifelse(
    !is.na(target_indices),
    paste0(value_list, ".", ave(seq_along(value_list), target_indices, FUN = seq_along)),
    value_list
  )
}

# Keep original wilcox testing method
perform_tests <- function(matrix1, matrix2) {
  results <- list()
  for (i in 1:nrow(matrix1)) {
    gene_id <- rownames(matrix1)[i]
    gene_result <- wilcox.test(as.numeric(matrix1[i, ]), as.numeric(matrix2[i, ]), paired = TRUE)
    gene_result$test_statistic <- gene_result$statistic
    results[[gene_id]] <- gene_result
  }
  return(results)
}

# ============================ #
#  LOOP THROUGH COMPARISONS   #
# ============================ #

comparison_list <- list(
  list(a = matrixhg18, b = matrixhg19, label = "GRCh37/NCBI36"),
  list(a = matrixhg18, b = matrixv30,  label = "GRCh38.p12/NCBI36"),
  list(a = matrixhg18, b = matrixv43,  label = "GRCh38.p13/NCBI36"),
  list(a = matrixhg18, b = matrixT2T,  label = "T2T-CHM13v2.0/NCBI36"),
  list(a = matrixhg19, b = matrixv30,  label = "GRCh38.p12/GRCh37"),
  list(a = matrixhg19, b = matrixv43,  label = "GRCh38.p13/GRCh37"),
  list(a = matrixhg19, b = matrixT2T,  label = "CHM13v2.0/GRCh37"),
  list(a = matrixv30,  b = matrixv43,  label = "GRCh38.p13/GRCh38.p12"),
  list(a = matrixv30,  b = matrixT2T,  label = "CHM13v2.0/GRCh38.p12"),
  list(a = matrixv43,  b = matrixT2T,  label = "CHM13v2.0/GRCh38.p13")
)

all_results <- list()

for (comp in comparison_list) {
  print(paste0("Running Wilcoxon test for: ", comp$label))
  a_common <- comp$a[rownames(comp$a) %in% rownames(comp$b), ]
  b_common <- comp$b[rownames(comp$b) %in% rownames(comp$a), ]

  gene_tests <- perform_tests(a_common, b_common)
  p_values <- sapply(gene_tests, function(x) x$p.value)
  v_values <- sapply(gene_tests, function(x) x$statistic)
  adj_p <- p.adjust(p_values, method = "bonferroni")

  results_df <- data.frame(
    Gene = names(gene_tests),
    P_Value = p_values,
    V_Value = v_values,
    adj.p.value = adj_p,
    Significant = ifelse(adj_p <= 0.05, "Yes", "No"),
    assembly = comp$label,
    stringsAsFactors = FALSE
  )
  all_results[[comp$label]] <- results_df
}

all2 <- do.call(rbind, all_results)
all2 <- all2[complete.cases(all2), ]

# ============================ #
#         PLOT RESULTS         #
# ============================ #

g <- ggplot(all2, aes(x = factor(assembly, levels = unique(assembly)), fill = Significant)) +
  geom_bar(stat = "count", position = "dodge") +
  labs(x = "Reference Comparisons", y = "Expressed Gene Count") +
  scale_fill_manual(
    values = c("red", "blue"),
    labels = c("> 0.05", "≤ 0.05"),
    name = "Bonferroni Correction"
  ) +
  theme_bw(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = c(0.05, 0.95),       # Top left corner
    legend.justification = c(0, 1),        # Align top-left of legend box
    legend.background = element_rect(fill = "white", color = "black") # Optional: makes it easier to see
  )

ggsave(
  filename = plot_path,
  plot = g,
  device = cairo_pdf,
  width = 10,
  height = 9
)

