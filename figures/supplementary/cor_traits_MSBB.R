# Supplementary Figure 1
# ==== Load Libraries ====
library(variancePartition)
library(CCP)

# ==== Load and Filter Data ====
infoall <- readRDS(file.path(data_dir, "info_all_v30_2023-10-23.RDS"))

infoall2 <- infoall[, sapply(infoall, function(col) !is.factor(col) || nlevels(col) >= 2)]

columns_to_remove <- c(
  "rest_id", "metabolomics.bileacids_id", "metabolomics.p180brain_id", 
  "metabolomics.p180serum_id", "age_at_visit", "age_first_ad_dx", 
  "cts_mmse30", "cts_mmse30_first_ad_dx", "dcfdx", "ad_reagan"
)
df <- infoall2[, !(names(infoall2) %in% columns_to_remove)]

# ==== Subset Traits and Define Formula ====
colnames(infoall2)[19]="BBscore"
subtraits <- infoall2[, c("CDR", "CERJ", "PlaqueMean", "BBscore")]
form <- ~ CDR + CERJ + PlaqueMean + BBscore

# ==== Canonical Correlation Analysis with Permutation Testing ====
canCorPairsTest <- function (formula, data, showWarnings = TRUE, number_perm = 1000) {
    formula <- stats::as.formula(formula)
    if (variancePartition:::.isMixedModelFormula(formula)) {
        stop("Invalid formula: mixed effects are not supported.")
    }

    X <- model.matrix(formula, data)
    varLabels <- attr(terms(formula), "term.labels")
    if (length(varLabels) < 2) {
        stop("Formula must include at least two variables.")
    }

    idx <- attr(X, "assign")
    variableList <- list()
    for (i in unique(idx)) {
        if (i == 0) next
        variableList[[varLabels[i]]] <- X[, idx == i, drop = FALSE]
    }

    checkNames <- names(which(sapply(variableList, sd) == 0))
    if (length(checkNames) != 0) {
        stop(paste("Variables with zero variance:", paste(checkNames, collapse = ", ")))
    }

    C <- matrix(NA, length(varLabels), length(varLabels), dimnames = list(varLabels, varLabels))
    diag(C) <- 1
    C_p <- C
    colinear_set <- c()
    pairs <- combn(varLabels, 2)

    for (i in 1:ncol(pairs)) {
        key1 <- pairs[1, i]
        key2 <- pairs[2, i]

        fit <- cancor(variableList[[key1]], variableList[[key2]])
        corr_mean <- mean(fit$cor)
        C[key1, key2] <- C[key2, key1] <- corr_mean

        capture.output({
            pv <- p.perm(variableList[[key1]], variableList[[key2]], nboot = number_perm)$p
        }, file = "/dev/null")

        C_p[key1, key2] <- C_p[key2, key1] <- pv

        if (showWarnings && corr_mean > 0.999) {
            colinear_set <- c(colinear_set, paste(key1, "and", key2))
        }
    }

    if (length(colinear_set) > 0) {
        warning("High colinearity detected:\n", paste(colinear_set, collapse = "\n"))
    }

    return(list(C = C, pvalues = C_p))
}

# ==== Run Analysis ====
number_perm <- 10000
ccFit <- canCorPairsTest(formula = form, data = subtraits, number_perm = number_perm)

# ==== P-value Adjustment and Plot Annotations ====
plotCorrMatrix(ccFit$C)

ccFit$pvalues_nonzero <- ccFit$pvalues
ccFit$pvalues_nonzero[ccFit$pvalues_nonzero == 0] <- 1 / (number_perm + 1)

adj.pvalues <- p.adjust(as.matrix(ccFit$pvalues_nonzero), method = "fdr")
dim(adj.pvalues) <- dim(ccFit$pvalues)
dimnames(adj.pvalues) <- dimnames(ccFit$pvalues)
ccFit$adj.pvalues <- adj.pvalues

ccFit$adj.pvalues_formatted <- formatC(ccFit$adj.pvalues, format = "e", digits = 1)
ccFit$adj.pvalues_formatted[ccFit$adj.pvalues > 0.05] <- ""

note <- paste(round(ccFit$C, 2), "\n(", ccFit$adj.pvalues_formatted, ")", sep = "")
note[ccFit$adj.pvalues_formatted == ""] <- round(ccFit$C[ccFit$adj.pvalues_formatted == ""], 2)
note[ccFit$pvalues == 0] <- paste(round(ccFit$C[ccFit$pvalues == 0], 2), "\n(<", ccFit$adj.pvalues_formatted[ccFit$pvalues == 0], ")", sep = "")
dim(note) <- dim(ccFit$pvalues)
ccFit$note <- note

# ==== Save Plot ====
pdf(file.path(plot_dir, "cor_traits_msbb.pdf"))
plotCorrMatrix(ccFit$C, cellnote = round(ccFit$C, 2), notecex = 1.2, notecol = "black")
dev.off()
