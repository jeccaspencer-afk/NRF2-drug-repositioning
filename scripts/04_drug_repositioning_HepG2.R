# --------------------------------------------------
# NRF2 Drug Repositioning Pipeline
# Author: Jessica Spencer
# Script: 04_drug_repositioning_HepG2.R
# Purpose: Score LINCS compound-induced signatures against the NRF2 consensus signature
# Inputs:
#   - results/NRF2_core14_consensus_TCGA_LIHC__LINCS_HEPG2.tsv
#   - data/LINCS/GSE92742_Broad_LINCS_sig_info.txt
#   - data/LINCS/GSE92742_Broad_LINCS_gene_info.txt
#   - data/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
# Outputs:
#   - paired_pr_gene_ids.txt
#   - compound_sig_scores_PKLRstyle.csv
#   - ranked_drugs_bestdose_PKLRstyle_inhibitors.csv
# --------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(cmapR)
  library(rhdf5)
})

cell <- "HEPG2"
chunk_size <- 1500
rank_for <- "inhibitors"

consensus_file <- "results/NRF2_core14_consensus_TCGA_LIHC__LINCS_HEPG2.tsv"
siginfo_txt <- "data/LINCS/GSE92742_Broad_LINCS_sig_info.txt"
geneinfo_txt <- "data/LINCS/GSE92742_Broad_LINCS_gene_info.txt"
level5_gctx <- "data/LINCS/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"

out_dir <- file.path(getwd(), paste0("NRF2_core14_finalstyle_", cell))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cons <- fread(consensus_file)
sig_info <- fread(siginfo_txt)
gene_info <- fread(geneinfo_txt)

if (!all(c("gene", "Z_combined") %in% names(cons))) {
  stop("Consensus file must contain columns: gene, Z_combined")
}

if (!all(c("pr_gene_id", "pr_gene_symbol") %in% names(gene_info))) {
  stop("gene_info must contain pr_gene_id and pr_gene_symbol")
}

sym2id <- setNames(as.character(gene_info$pr_gene_id), gene_info$pr_gene_symbol)

# Prepare consensus vector
cons[, pr_gene_id := unname(sym2id[gene])]
cons <- cons[!is.na(pr_gene_id) & pr_gene_id != ""]

cons[, absZ := abs(Z_combined)]
setorder(cons, pr_gene_id, -absZ)
cons <- cons[!duplicated(pr_gene_id)]
cons[, absZ := NULL]

gctx_rids <- h5read(level5_gctx, "/0/META/ROW/id")
keep_ids <- intersect(cons$pr_gene_id, gctx_rids)

cons <- cons[pr_gene_id %in% keep_ids]
cat("Consensus genes after mapping and intersection:", nrow(cons), "\n")

if (nrow(cons) < 50) {
  stop("Too few overlapping genes between consensus signature and GCTX matrix.")
}

cons_vec <- cons$Z_combined
names(cons_vec) <- cons$pr_gene_id

# Select compound signatures for chosen cell line
req_sig <- c("sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_dose", "pert_time")
miss <- setdiff(req_sig, names(sig_info))
if (length(miss) > 0) {
  stop(paste("sig_info missing:", paste(miss, collapse = ", ")))
}

cp_meta <- sig_info[pert_type == "trt_cp" & cell_id == cell]
cp_ids <- unique(as.character(cp_meta$sig_id))
cp_ids <- cp_ids[!is.na(cp_ids) & cp_ids != ""]

cat("Compound signatures in ", cell, ": ", length(cp_ids), "\n", sep = "")
if (length(cp_ids) == 0) {
  stop("No compound signatures found for this cell line.")
}

# Spearman correlation in chunks
rid_use <- names(cons_vec)
n_genes <- length(rid_use)

cat("Genes used in correlation:", n_genes, "\n")
writeLines(rid_use, file.path(out_dir, "paired_pr_gene_ids.txt"))

n_chunks <- ceiling(length(cp_ids) / chunk_size)
cat("Scoring in ", n_chunks, " chunks\n", sep = "")

all_scores <- vector("list", n_chunks)

for (i in seq_len(n_chunks)) {
  idx <- ((i - 1) * chunk_size + 1):min(i * chunk_size, length(cp_ids))
  chunk_cids <- cp_ids[idx]
  
  cat("Chunk ", i, "/", n_chunks, " (", length(chunk_cids), " signatures)\n", sep = "")
  mat <- parse_gctx(level5_gctx, rid = rid_use, cid = chunk_cids)@mat
  
  r <- as.numeric(cor(cons_vec, mat, method = "spearman", use = "pairwise.complete.obs"))
  names(r) <- colnames(mat)
  
  tstat <- r * sqrt((n_genes - 2) / pmax(1 - r^2, 1e-16))
  pval <- 2 * pt(-abs(tstat), df = n_genes - 2)
  
  all_scores[[i]] <- data.table(sig_id = names(r), r = r, p = pval)
  rm(mat)
  gc()
}

sig_scores <- rbindlist(all_scores)
sig_scores[, FDR := p.adjust(p, method = "BH")]

fwrite(sig_scores, file.path(out_dir, "compound_sig_scores_PKLRstyle.csv"))

# Collapse to best dose/time per drug
# For inhibitors, retain the most negative correlation per compound.
# For activators, retain the most positive correlation per compound.
cp_meta2 <- cp_meta[, .(sig_id, pert_id, pert_iname, pert_dose, pert_time)]
m <- merge(sig_scores, cp_meta2, by = "sig_id")

if (rank_for == "inhibitors") {
  setorder(m, pert_id, r)
  best <- m[, .SD[1], by = pert_id]
  setorder(best, r)
} else {
  setorder(m, pert_id, -r)
  best <- m[, .SD[1], by = pert_id]
  setorder(best, -r)
}

best <- merge(best, m[, .(n_sigs = .N), by = pert_id], by = "pert_id")

if (rank_for == "inhibitors") {
  setorder(best, r)
} else {
  setorder(best, -r)
}

fwrite(best, file.path(out_dir, paste0("ranked_drugs_bestdose_PKLRstyle_", rank_for, ".csv")))

cat("Unique compounds ranked:", nrow(best), "\n")
cat("Compounds with FDR < 0.05:", sum(best$FDR < 0.05, na.rm = TRUE), "\n")

cat("\nTop 20 ", rank_for, " (Spearman connectivity):\n", sep = "")
print(best[1:20, .(pert_id, pert_iname, r, FDR, n_sigs, pert_dose, pert_time)])

# Sanity check
check_terms <- c("tert-butylhydroquinone", "tBHQ", "sulforaphane", "bardoxolone")
chk <- best[Reduce(`|`, lapply(check_terms, function(x) grepl(x, pert_iname, ignore.case = TRUE)))]

if (nrow(chk) > 0) {
  cat("\nSanity check hits:\n")
  print(chk[, .(pert_iname, r, FDR, pert_dose, pert_time, n_sigs)])
}

cat("\nDONE. Outputs written to:\n", out_dir, "\n")
