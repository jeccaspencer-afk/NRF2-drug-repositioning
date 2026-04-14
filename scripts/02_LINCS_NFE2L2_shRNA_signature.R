# --------------------------------------------------
# NRF2 Drug Repositioning Pipeline
# Author: Jessica Spencer
# Script: 02_LINCS_NFE2L2_shRNA_signature.R
# Purpose: Construct LINCS-derived NFE2L2 knockdown signature
# Inputs:
#   - GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
#   - GSE92742_Broad_LINCS_sig_info.txt
#   - GSE92742_Broad_LINCS_gene_info.txt
# Outputs:
#   - NRF2_LINCS_NFE2L2_shRNA_consensus_full_*.tsv
#   - NRF2_LINCS_NFE2L2_KD_DOWN_top*.tsv
#   - NRF2_LINCS_NFE2L2_KD_UP_top*.tsv
# --------------------------------------------------

library(data.table)
library(dplyr)
library(cmapR)
library(matrixStats)

gctx_file <- "data/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
sig_info_file <- "data/LINCS/GSE92742_Broad_LINCS_sig_info.txt"
gene_info_file <- "data/LINCS/GSE92742_Broad_LINCS_gene_info.txt"

dir.create("results", showWarnings = FALSE)

USE_HEPG2_ONLY <- TRUE
top_n <- 200

# Load metadata
sig_info <- fread(sig_info_file, data.table = FALSE)
gene_info <- fread(gene_info_file, data.table = FALSE)

# Select NFE2L2 shRNA signatures
nrf2_sh <- sig_info %>%
  filter(pert_type == "trt_sh", pert_iname == "NFE2L2")

if (USE_HEPG2_ONLY) {
  nrf2_sh <- nrf2_sh %>% filter(cell_id == "HEPG2")
}

cat("NFE2L2 trt_sh signatures:", nrow(nrf2_sh), "\n")
cat("Cell lines represented:\n")
print(sort(table(nrf2_sh$cell_id), decreasing = TRUE))

if (nrow(nrf2_sh) == 0) {
  stop("No NFE2L2 shRNA signatures found with current filters.")
}

# Load selected signatures
mat <- parse_gctx(gctx_file, cid = nrf2_sh$sig_id)
z <- mat@mat

# Map probe IDs to gene symbols
map <- unique(gene_info[, c("pr_gene_id", "pr_gene_symbol")])
row_syms <- map$pr_gene_symbol[match(mat@rid, map$pr_gene_id)]
row_syms[is.na(row_syms)] <- mat@rid[is.na(row_syms)]

# Compute median knockdown signature
cons_z <- rowMedians(z, na.rm = TRUE)

cons_df <- data.frame(
  gene = row_syms,
  z_median = cons_z,
  stringsAsFactors = FALSE
) %>%
  group_by(gene) %>%
  summarise(z_median = median(z_median, na.rm = TRUE), .groups = "drop") %>%
  arrange(z_median)

# Create top up/down gene sets
NRF2_KD_DOWN <- cons_df %>% slice_head(n = top_n)
NRF2_KD_UP <- cons_df %>% arrange(desc(z_median)) %>% slice_head(n = top_n)

# Save outputs
tag <- if (USE_HEPG2_ONLY) "HEPG2" else "ALLCELLS"

write.table(
  cons_df,
  paste0("results/NRF2_LINCS_NFE2L2_shRNA_consensus_full_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  NRF2_KD_DOWN,
  paste0("results/NRF2_LINCS_NFE2L2_KD_DOWN_top", top_n, "_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  NRF2_KD_UP,
  paste0("results/NRF2_LINCS_NFE2L2_KD_UP_top", top_n, "_", tag, ".tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Sanity check
targets <- c(
  "ABHD4", "AKR1C3", "EPHX1", "FTH1", "FTL", "GCLC", "GCLM",
  "GSR", "ME1", "NQO1", "OSGIN1", "PIR", "SLC7A11", "SRXN1"
)

cat("\nCurated NRF2 target genes in LINCS NFE2L2 knockdown signature:\n")
print(cons_df %>% filter(gene %in% targets) %>% arrange(z_median))
