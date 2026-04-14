# --------------------------------------------------
# NRF2 Drug Repositioning Pipeline
# Author: Jessica Spencer
# Script: 01_TCGA_NRF2_coexpression_signature.R
# Purpose: Construct TCGA-LIHC NRF2 co-expression signature
# Inputs: HPA_V2.TCGA.LIHC.data.tpm.tsv
# Outputs:
#   - NRF2_LIHC_positive_*.tsv
#   - NRF2_LIHC_negative_*.tsv
#   - NRF2_LIHC_full_ranked_*.tsv
# --------------------------------------------------

library(data.table)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)

tpm_file <- "data/HPA_V2.TCGA.LIHC.data.tpm.tsv"

min_tpm <- 1
min_frac_samples <- 0.20
fdr_cutoff <- 0.01
rho_cutoff <- 0.30

p_floor <- 1e-300
p_cap <- 1 - 1e-16

targets <- c(
  "ABHD4", "AKR1C3", "EPHX1", "FTH1", "FTL", "GCLC", "GCLM",
  "GSR", "ME1", "NQO1", "OSGIN1", "PIR", "SLC7A11", "SRXN1"
)

dir.create("results", showWarnings = FALSE)

# Load TCGA-LIHC TPM
tpm <- fread(tpm_file, data.table = FALSE)

expr_cols <- grepl("^ENSG", colnames(tpm))
expr <- tpm[, expr_cols]
meta <- tpm[, !expr_cols]

stopifnot(ncol(expr) > 1000)
stopifnot(nrow(expr) > 100)

# Filter low-expression genes
keep_genes <- colSums(expr >= min_tpm) >= (min_frac_samples * nrow(expr))
expr_filt <- expr[, keep_genes]

# Log transform
expr_log <- log2(expr_filt + 1)

# Map ENSG to SYMBOL
ensg_ids <- colnames(expr_log)
ensg_clean <- sub("\\..*$", "", ensg_ids)

symbol_map <- suppressMessages(
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(ensg_clean),
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  )
)

symbol_map <- symbol_map[!is.na(symbol_map$SYMBOL), ]
symbol_map <- symbol_map[!duplicated(symbol_map$ENSEMBL), ]

ensg2sym <- setNames(symbol_map$SYMBOL, symbol_map$ENSEMBL)
gene_symbols <- ensg2sym[ensg_clean]

# Compute NRF2 activity score
sym2idx <- split(seq_along(gene_symbols), gene_symbols)
target_idx <- unlist(sym2idx[targets], use.names = FALSE)
target_idx <- target_idx[!is.na(target_idx)]

if (length(target_idx) < 3) {
  stop("Too few NRF2 target genes found after ENSG->SYMBOL mapping.")
}

nrf2_score <- rowMeans(expr_log[, target_idx, drop = FALSE], na.rm = TRUE)

# Spearman correlation
rho <- cor(expr_log, nrf2_score, method = "spearman")

n <- nrow(expr_log)
t_stat <- rho * sqrt((n - 2) / pmax(1 - rho^2, 1e-12))
pval <- 2 * pt(-abs(t_stat), df = n - 2)
pval <- pmin(pmax(pval, p_floor), p_cap)

cor_results <- data.frame(
  ensg = colnames(expr_log),
  symbol = gene_symbols,
  rho = as.numeric(rho),
  p = as.numeric(pval),
  stringsAsFactors = FALSE
)

# FDR and signed Z-score
cor_results$FDR <- p.adjust(cor_results$p, method = "BH")

p_use <- pmin(pmax(cor_results$p, p_floor), p_cap)
cor_results$Z_signed <- qnorm(p_use / 2, lower.tail = FALSE) * sign(cor_results$rho)

# Strict NRF2 signatures
NRF2_positive <- cor_results %>%
  filter(rho >= rho_cutoff, FDR < fdr_cutoff) %>%
  arrange(desc(Z_signed))

NRF2_negative <- cor_results %>%
  filter(rho <= -rho_cutoff, FDR < fdr_cutoff) %>%
  arrange(Z_signed)

# Save outputs
tag <- paste0("FDR", fdr_cutoff, "_rho", rho_cutoff)

write.table(
  NRF2_positive,
  paste0("results/NRF2_LIHC_positive_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  NRF2_negative,
  paste0("results/NRF2_LIHC_negative_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  cor_results,
  paste0("results/NRF2_LIHC_full_ranked_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Summary
cat("NRF2 targets found (symbols):\n")
print(intersect(targets, unique(cor_results$symbol)))

cat("\nNumber of curated NRF2 target genes found:",
    length(intersect(targets, unique(cor_results$symbol))),
    "of", length(targets), "\n")

cat("\nCounts:\n")
cat("  Positive genes:", nrow(NRF2_positive), "\n")
cat("  Negative genes:", nrow(NRF2_negative), "\n")

cat("\nAny infinite Z-scores?\n")
print(sum(is.infinite(cor_results$Z_signed)))

cat("\nCanonical targets that pass strict thresholds:\n")
print(NRF2_positive %>% filter(symbol %in% targets) %>% arrange(desc(Z_signed)))
