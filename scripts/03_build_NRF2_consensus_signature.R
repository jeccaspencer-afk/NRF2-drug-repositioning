# --------------------------------------------------
# NRF2 Drug Repositioning Pipeline
# Author: Jessica Spencer
# Script: 03_build_NRF2_consensus_signature.R
# Purpose: Integrate TCGA co-expression and LINCS NFE2L2 knockdown signatures
# Inputs:
#   - NRF2_LIHC_full_ranked_FDR0.01_rho0.3.tsv
#   - NRF2_LINCS_NFE2L2_shRNA_consensus_full_HEPG2.tsv
# Outputs:
#   - NRF2_consensus_TCGA_LIHC_LINCS_HEPG2.tsv
#   - NRF2_consensus_UP_top200.tsv
#   - NRF2_consensus_DOWN_top200.tsv
# --------------------------------------------------

library(data.table)

tcga_file  <- "results/NRF2_LIHC_full_ranked_FDR0.01_rho0.3.tsv"
lincs_file <- "results/NRF2_LINCS_NFE2L2_shRNA_consensus_full_HEPG2.tsv"

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# Load data
tcga <- fread(tcga_file)
lincs <- fread(lincs_file)

# TCGA layer
if (all(c("symbol", "p", "rho") %in% colnames(tcga))) {
  p_use <- pmin(pmax(tcga$p, 1e-300), 1 - 1e-16)
  
  tcga_layer <- data.table(
    gene   = tcga$symbol,
    Z_TCGA = qnorm(p_use / 2, lower.tail = FALSE) * sign(tcga$rho)
  )
  
} else if (all(c("symbol", "Z_signed") %in% colnames(tcga))) {
  tcga_layer <- data.table(
    gene   = tcga$symbol,
    Z_TCGA = tcga$Z_signed
  )
  
  tcga_layer[is.infinite(Z_TCGA), Z_TCGA := NA_real_]
  
  z_cap <- 12
  tcga_layer[, Z_TCGA := pmax(pmin(Z_TCGA, z_cap), -z_cap)]
  
} else {
  stop("TCGA file must contain either (symbol, p, rho) or (symbol, Z_signed).")
}

tcga_layer <- tcga_layer[!is.na(gene) & gene != "" & !is.na(Z_TCGA)]

tcga_layer[, absZ := abs(Z_TCGA)]
setorder(tcga_layer, gene, -absZ)
tcga_layer <- tcga_layer[!duplicated(gene)]
tcga_layer[, absZ := NULL]

# LINCS layer
stopifnot(all(c("gene", "z_median") %in% colnames(lincs)))

lincs_layer <- data.table(
  gene    = lincs$gene,
  Z_LINCS = -lincs$z_median
)

lincs_layer <- lincs_layer[!is.na(gene) & gene != "" & !is.na(Z_LINCS)]

lincs_layer[, absZ := abs(Z_LINCS)]
setorder(lincs_layer, gene, -absZ)
lincs_layer <- lincs_layer[!duplicated(gene)]
lincs_layer[, absZ := NULL]

# Merge layers
m <- merge(tcga_layer, lincs_layer, by = "gene")
cat("Genes in merged set:", nrow(m), "\n")

# Directional concordance
consensus <- m[
  (Z_TCGA > 0 & Z_LINCS > 0) |
    (Z_TCGA < 0 & Z_LINCS < 0)
]

cat("Genes passing directional concordance:", nrow(consensus), "\n")

# Combined Z-score
consensus[, Z_combined := (Z_TCGA + Z_LINCS) / sqrt(2)]
setorder(consensus, -Z_combined)

cat("Any infinite combined Z? ", sum(is.infinite(consensus$Z_combined)), "\n")

cat("Consensus genes with positive combined Z:", sum(consensus$Z_combined > 0), "\n")
cat("Consensus genes with negative combined Z:", sum(consensus$Z_combined < 0), "\n")

# Save outputs
f_full <- file.path(out_dir, "NRF2_consensus_TCGA_LIHC__LINCS_HEPG2.tsv")
fwrite(consensus, f_full, sep = "\t")

top_n <- 200

NRF2_UP <- consensus[Z_combined > 0][1:min(top_n, sum(consensus$Z_combined > 0))]

NRF2_DOWN <- consensus[Z_combined < 0]
setorder(NRF2_DOWN, Z_combined)
NRF2_DOWN <- NRF2_DOWN[1:min(top_n, nrow(NRF2_DOWN))]

f_up <- file.path(out_dir, paste0("NRF2_consensus_UP_top", top_n, ".tsv"))
f_dn <- file.path(out_dir, paste0("NRF2_consensus_DOWN_top", top_n, ".tsv"))

fwrite(NRF2_UP, f_up, sep = "\t")
fwrite(NRF2_DOWN, f_dn, sep = "\t")

targets <- c(
  "ABHD4", "AKR1C3", "EPHX1", "FTH1", "FTL", "GCLC", "GCLM",
  "GSR", "ME1", "NQO1", "OSGIN1", "PIR", "SLC7A11", "SRXN1"
)

cat("\nCurated NRF2 target genes in consensus:\n")
print(
  consensus[gene %in% targets, .(gene, Z_TCGA, Z_LINCS, Z_combined)][order(-Z_combined)]
)

cat("\nNumber of curated NRF2 target genes found in consensus:",
    sum(targets %in% consensus$gene), "of", length(targets), "\n")
