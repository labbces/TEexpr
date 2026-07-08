#!/usr/bin/env Rscript
# ==============================================================================
# pearson_coexpression_network.R
#
# PURPOSE:
#   Build a gene/TE coexpression network by computing a full pairwise Pearson
#   correlation matrix from RUVr-normalized count data.
#
#   Computation is done in ROW CHUNKS, processed in parallel (mclapply) and
#   written to disk incrementally, instead of holding the full n_genes x
#   n_genes correlation matrix in memory at once. This mirrors the strategy
#   used in the sugarcane pipeline script for large sample/gene counts, and
#   keeps memory use controllable regardless of how many samples a given
#   species has.
#
# ASSUMPTIONS (important):
#   - Input is one of the OUTPUTS of ruv_analysis.R:
#       RUVr_k<k>_<species>_corrected_counts.tsv   (RUVr-corrected counts), or
#       RUVr_k<k>_<species>_vst.tsv                (VST of the RUVr-corrected counts)
#     CV filtering (CV >= 15%, cv_threshold in ruv_analysis.R) is applied to
#     RAW counts BEFORE RUVr correction in that script (see filter_by_cv()),
#     so genes/TEs here are already restricted to the CV-filtered set.
#     NO additional CV filtering is applied here, to avoid recomputing a
#     decision already made upstream.
#   - Note: CV was computed on raw counts, not on VST. If you use the VST
#     file as input here (recommended for Pearson correlation, since VST
#     stabilizes variance across expression levels), keep in mind the gene
#     set itself was decided on the raw-count scale.
#   - If you later build a network on a SUBSET of samples (e.g. only one
#     species, one tissue, one condition) that differs from the sample set
#     used to compute CV in ruv_analysis.R, you should re-evaluate whether
#     CV filtering needs to be recalculated for that subset.
#
# INPUT:
#   A tab-separated matrix of RUVr-normalized counts, genes/TEs in rows,
#   samples in columns, with row names = GeneID/TE_ID and header = sample IDs.
#   Example: RUVr_k1_Sbicolor_vst.tsv
#
# OUTPUTS:
#   1. matrix_<label>_pearson.tsv    - full pairwise Pearson correlation matrix
#   2. matrix_<label>_pvalues.tsv    - full pairwise RAW p-value matrix (t-test on r)
#   3. edges_<label>_r<r>_q<q>.tsv   - filtered edge list (|r| > r AND FDR q < q),
#                                      ready for Cytoscape / igraph import
#
# NOTE ON MULTIPLE TESTING (FDR / Benjamini-Hochberg):
#   With n_genes/TEs, the number of pairwise tests is n*(n-1)/2 — often in
#   the tens of millions. Using a raw p-value cutoff (e.g. p < 0.05) without
#   correction means ~5% of ALL tested pairs are expected to look
#   "significant" purely by chance, even if no real association exists
#   anywhere (e.g. 50 million pairs -> ~2.5 million false positives from
#   noise alone). The Benjamini-Hochberg (BH) procedure corrects for this by
#   adjusting p-values based on their rank among ALL tests performed, giving
#   a q-value: the expected proportion of false positives among everything
#   called "significant" at that threshold. This requires seeing every raw
#   p-value from the WHOLE dataset at once (correction cannot be done
#   correctly per chunk in isolation, since it depends on the total number
#   of tests and their global ranking) — hence the dedicated step 6 below,
#   which runs only after every chunk has finished.
#
# NOTE ON SCOPE:
#   This script ONLY computes correlation (the coexpression "backbone").
#   Topological Overlap Matrix (TOM) and module detection (WGCNA) are
#   deliberately kept in a SEPARATE downstream script, so this correlation
#   matrix can be reused without recomputation (see chat discussion).
# ==============================================================================

library(data.table)   # fast I/O and table handling
library(matrixStats)  # rowSds(), used only for a quick sanity-check summary
library(parallel)     # mclapply() for chunk-level parallelism

# ==============================================================================
# COMMAND-LINE ARGUMENTS
# Cluster: Rscript pearson_coexpression_network.R --in_file /path --label Sbicolor ...
# Interactive: fallback values in CONFIGURATION section are used
# (same parse_arg() convention as ruv_analysis.R, for consistency)
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else NULL
}

cli_in_file    <- parse_arg(args, "--in_file")
cli_out_dir    <- parse_arg(args, "--out_dir")
cli_label      <- parse_arg(args, "--label")
cli_chunk_rows <- parse_arg(args, "--chunk_rows")
cli_n_cores    <- parse_arg(args, "--n_cores")
cli_r_thresh   <- parse_arg(args, "--r_thresh")
cli_q_thresh   <- parse_arg(args, "--q_thresh")
cli_te_pattern <- parse_arg(args, "--te_pattern")

# ==============================================================================
# CONFIGURATION - CLI values take priority; fallbacks used for interactive runs
# ==============================================================================
in_file    <- if (!is.null(cli_in_file))    cli_in_file    else "~/pearson_test/RUVr_k1_Sitalica_vst.tsv"
out_dir    <- if (!is.null(cli_out_dir))    cli_out_dir    else "~/pearson_test/italica_results"
label      <- if (!is.null(cli_label))      cli_label      else "italica"
CHUNK_ROWS <- if (!is.null(cli_chunk_rows)) as.integer(cli_chunk_rows) else 500
N_CORES    <- if (!is.null(cli_n_cores))    as.integer(cli_n_cores)    else 50
R_THRESH   <- if (!is.null(cli_r_thresh))   as.numeric(cli_r_thresh)   else 0.8
Q_THRESH   <- if (!is.null(cli_q_thresh))   as.numeric(cli_q_thresh)   else 0.05  # FDR-adjusted (BH) threshold

# Regex to flag TE vs gene IDs in the edge list metadata.
# Matches the tx2gene convention used in ruv_analysis.R, where all TE
# instances/representative elements are prefixed with "TE_" and genes are not.
te_id_pattern <- if (!is.null(cli_te_pattern)) cli_te_pattern else "^TE_"

# mclapply() parallelizes via fork(), which is a Unix/Linux mechanism and is
# NOT available on Windows. This only matters for local testing in RStudio
# on a Windows machine — on the SGE cluster (Linux) this is a non-issue.
# We fall back to serial execution (mc.cores = 1) automatically so a local
# test run on Windows still works, just without real parallelism; use the
# actual cluster submission for real multi-core runs.
if (.Platform$OS.type == "windows" && N_CORES > 1) {
  warning("mclapply() does not support mc.cores > 1 on Windows. ",
          "Falling back to mc.cores = 1 for this local test run. ",
          "Submit via the SGE cluster script for real parallel execution.")
  N_CORES <- 1L
}

# Prevent BLAS from spawning its own internal threads on top of the mclapply
# workers below — without this, N_CORES workers x BLAS-internal-threads can
# massively oversubscribe the machine's CPUs and slow everything down.
Sys.setenv(OMP_NUM_THREADS      = 1,
           OPENBLAS_NUM_THREADS = 1,
           MKL_NUM_THREADS      = 1,
           BLAS_NUM_THREADS     = 1)

cat("================================================================\n")
cat(sprintf("PEARSON COEXPRESSION NETWORK — label: %s\n", label))
cat(sprintf("Input file: %s\n", in_file))
cat(sprintf("Output dir: %s\n", out_dir))
cat(sprintf("Chunk rows: %d | Workers: %d | r threshold: %.2f | FDR q threshold: %.3f\n",
            CHUNK_ROWS, N_CORES, R_THRESH, Q_THRESH))
cat("================================================================\n")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
message("Loading RUVr-normalized count matrix: ", in_file)
vst <- as.matrix(read.table(in_file, sep = "\t", header = TRUE, row.names = 1))

message("Matrix dimensions (genes/TEs x samples): ",
        paste(dim(vst), collapse = " x "))

# Quick sanity check: any NA or zero-variance rows would break correlation
n_na <- sum(is.na(vst))
if (n_na > 0) {
  warning(n_na, " NA values found in the input matrix. ",
          "Check upstream RUVr output before proceeding.")
}

row_sd <- rowSds(vst)
n_zero_var <- sum(row_sd == 0, na.rm = TRUE)
if (n_zero_var > 0) {
  warning(n_zero_var, " genes/TEs have zero variance across samples ",
          "(constant values) and will produce NA correlations. ",
          "Consider removing them before proceeding.")
}

genes   <- rownames(vst)
n_genes <- length(genes)
n_samp  <- ncol(vst)
vst_t   <- t(vst)  # samples x genes, reused as the "against all genes" side of each block
df      <- n_samp - 2L  # degrees of freedom for the correlation t-test below

n_chunks    <- ceiling(n_genes / CHUNK_ROWS)
total_pairs <- n_genes * (n_genes - 1) / 2

message("Genes/TEs: ", format(n_genes, big.mark = ","))
message("Samples: ", n_samp)
message("Chunks: ", n_chunks, " | Workers: ", N_CORES)
message("Total pairwise tests (upper triangle): ", format(total_pairs, big.mark = ","))
message(sprintf("Approx. RAM per in-flight block: ~%.2f GB (%d workers x %d rows x %s genes x 8 bytes)",
                (N_CORES * CHUNK_ROWS * n_genes * 8) / 1e9, N_CORES, CHUNK_ROWS,
                format(n_genes, big.mark = ",")))

# ==============================================================================
# 2. TEMP DIRECTORY FOR CHUNK OUTPUTS
# ==============================================================================
tmp_dir <- file.path(out_dir, paste0(".tmp_", label))
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 3. WORKER: one chunk = CHUNK_ROWS genes/TEs correlated against ALL genes/TEs
# ==============================================================================
# Each worker only ever holds a (CHUNK_ROWS x n_genes) block in memory, never
# the full (n_genes x n_genes) matrix — this is what keeps RAM controllable.
#
# For each correlation value, a RAW p-value is derived via the standard
# t-test for Pearson correlation: t = r * sqrt(df / (1 - r^2)), with
# df = n_samp - 2. This tests H0: rho = 0 (no linear association).
#
# IMPORTANT: workers no longer decide which pairs are "significant" on their
# own. Correct FDR correction (see note at the top of this script) needs the
# full set of raw p-values across ALL pairs, so each worker writes out every
# upper-triangle pair from its block (Node1, Node2, r, raw p) UNFILTERED —
# not just the ones passing a threshold. Filtering by |r| and by the
# FDR-adjusted q-value happens later, in step 6, once every chunk is done.
worker_chunk <- function(chunk_idx, row_start, row_end) {
  tryCatch({
    chunk_genes <- genes[row_start:row_end]
    vst_chunk   <- vst[row_start:row_end, , drop = FALSE]
    
    cor_block <- cor(t(vst_chunk), vst_t)
    
    # --- derive p-values from correlation via the t-test above -------------
    # 1e-15 avoids division by zero for pairs with |r| == 1 (e.g. self-pairs)
    t_stat     <- cor_block * sqrt(df / (1 - cor_block^2 + 1e-15))
    pval_block <- 2 * pt(-abs(t_stat), df = df)
    
    cor_block_r  <- round(cor_block, 3)
    pval_block_s <- signif(pval_block, 4)
    
    # --- write this block's rows to temp files (for the full matrices) ----
    tmp_cor  <- file.path(tmp_dir, sprintf("chunk_%05d_cor.tsv",  chunk_idx))
    tmp_pval <- file.path(tmp_dir, sprintf("chunk_%05d_pval.tsv", chunk_idx))
    
    cor_dt  <- data.table(GeneID = chunk_genes, as.data.table(cor_block_r))
    pval_dt <- data.table(GeneID = chunk_genes, as.data.table(pval_block_s))
    setnames(cor_dt,  c("GeneID", genes))
    setnames(pval_dt, c("GeneID", genes))
    
    fwrite(cor_dt,  file = tmp_cor,  sep = "\t", quote = FALSE, col.names = FALSE)
    fwrite(pval_dt, file = tmp_pval, sep = "\t", quote = FALSE, col.names = FALSE)
    
    # --- write this block's UNFILTERED upper-triangle pairs ----------------
    # (needed for the global FDR correction in step 6 — see note above)
    row_global_idx <- row_start:row_end
    # mask[i, j] = TRUE only if global column j is "after" global row i,
    # i.e. keeps each pair once (upper triangle) even across chunk boundaries
    mask <- outer(row_global_idx, seq_len(n_genes), FUN = function(r, c) c > r)
    
    idx <- which(mask, arr.ind = TRUE)
    
    tmp_pairs <- file.path(tmp_dir, sprintf("chunk_%05d_pairs.tsv", chunk_idx))
    pairs_dt <- data.table(
      Node1  = chunk_genes[idx[, 1]],
      Node2  = genes[idx[, 2]],
      weight = cor_block_r[idx],
      pvalue = pval_block_s[idx]
    )
    fwrite(pairs_dt, file = tmp_pairs, sep = "\t", quote = FALSE, col.names = FALSE)
    
    rm(cor_block, t_stat, pval_block, cor_block_r, pval_block_s,
       vst_chunk, cor_dt, pval_dt, mask, idx, pairs_dt); gc()
    
    list(ok = TRUE, cor = tmp_cor, pval = tmp_pval, pairs = tmp_pairs, chunk = chunk_idx)
    
  }, error = function(e) {
    list(ok = FALSE, chunk = chunk_idx, msg = conditionMessage(e))
  })
}

# ==============================================================================
# 4. CHUNK INDEX + PARALLEL COMPUTATION
# ==============================================================================
chunk_starts <- seq(1L, n_genes, by = CHUNK_ROWS)
chunk_ends   <- pmin(chunk_starts + CHUNK_ROWS - 1L, n_genes)
chunk_ids    <- seq_along(chunk_starts)

message("Launching parallel computation...")
results <- mclapply(
  chunk_ids,
  function(ci) worker_chunk(ci, chunk_starts[ci], chunk_ends[ci]),
  mc.cores           = N_CORES,
  mc.preschedule     = FALSE,
  mc.allow.recursive = FALSE
)

# --- check for errors before touching any output --------------------------
failed <- which(!vapply(results, `[[`, logical(1L), "ok"))
if (length(failed) > 0L) {
  msgs <- vapply(results[failed], `[[`, character(1L), "msg")
  stop("Errors in chunks ", paste(failed, collapse = ", "), ":\n",
       paste(msgs, collapse = "\n"))
}
message("All ", n_chunks, " chunks computed successfully.")

results <- results[order(vapply(results, `[[`, integer(1L), "chunk"))]

# ==============================================================================
# 5. CONCATENATE FULL CORRELATION AND RAW P-VALUE MATRICES
# ==============================================================================
out_matrix <- file.path(out_dir, paste0("matrix_", label, "_pearson.tsv"))
out_pvals  <- file.path(out_dir, paste0("matrix_", label, "_pvalues.tsv"))
header_line <- paste(c("GeneID", genes), collapse = "\t")
writeLines(header_line, con = out_matrix)
writeLines(header_line, con = out_pvals)

for (res in results) {
  file.append(out_matrix, res$cor)
  file.append(out_pvals,  res$pval)
}
message("Saved full correlation matrix: ", out_matrix)
message("Saved full RAW p-value matrix: ", out_pvals)

# ==============================================================================
# 6. GLOBAL FDR (BENJAMINI-HOCHBERG) CORRECTION + FILTERED EDGE LIST
# ==============================================================================
# This is the step that requires ALL pairs to be seen together (see the note
# at the top of the script). We read back every chunk's unfiltered pairs
# file, combine them into one table covering every tested pair exactly once,
# then apply p.adjust(method = "BH") across the WHOLE pvalue column at once.
# Only after that do we filter by |r| and by the FDR-adjusted q-value.
message("Loading all pairwise tests for global FDR correction (", 
        format(total_pairs, big.mark = ","), " pairs)...")

all_pairs <- rbindlist(
  lapply(results, function(res) {
    fread(res$pairs, header = FALSE,
          col.names = c("Node1", "Node2", "weight", "pvalue"))
  })
)

message("Applying Benjamini-Hochberg (BH) correction across all tests...")
all_pairs[, qvalue := p.adjust(pvalue, method = "BH")]

message("Filtering edges (|r| > ", R_THRESH, " and FDR q < ", Q_THRESH, ")...")
edge_list <- all_pairs[abs(weight) > R_THRESH & qvalue < Q_THRESH]
rm(all_pairs); gc()

# Flag each node as "gene" or "TE" based on the ID pattern above.
# This makes it easy downstream to separate TE-TE, TE-gene, and gene-gene
# edges for interpretation.
edge_list[, Node1_type := ifelse(grepl(te_id_pattern, Node1), "TE", "gene")]
edge_list[, Node2_type := ifelse(grepl(te_id_pattern, Node2), "TE", "gene")]

out_edges <- file.path(out_dir, paste0("edges_", label, "_r", R_THRESH, "_q", Q_THRESH, ".tsv"))
fwrite(edge_list, file = out_edges, sep = "\t", quote = FALSE)
message("Saved filtered edge list: ", out_edges)
message("Total edges above threshold: ", nrow(edge_list))

# ==============================================================================
# 7. CLEANUP
# ==============================================================================
unlink(tmp_dir, recursive = TRUE)

message("\nDone. Next step (separate script): soft-power selection + TOM + module detection.")
