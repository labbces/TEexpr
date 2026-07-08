# ==============================================================================
# Quantification Uncertainty Analysis Using Salmon Gibbs Replicates
#
# Compatible with single and multiple samples.
#
# Description:
#   Imports Salmon output (with Gibbs replicates) via tximport, builds a
#   SummarizedExperiment, scales inferential replicates, and computes five
#   uncertainty metrics — CV, SD, IQR, quasi-Poisson phi, and InfRV — across
#   all features (TEs and protein-coding genes) for all samples and conditions.
#   TE copies receive family/superfamily classification from the TEdistill
#   *.mapids file; genes receive feature_type = "gene".
#   Visualisation is handled separately by uncertainty_plots.R, which reads
#   the tables produced here.
#
# Input files:
#   CONDITIONS_FILE     -- TSV with columns: sample, condition
#   SALMON_DIR          -- Directory containing <sample>/quant.sf files
#   CLASSIFICATION_FILE -- TEdistill *.mapids file (no header, 3 columns:
#                          te_id | classification | te_id#classification)
#
# Output files (written to OUT_DIR/tables/):
#   01_run_summary.tsv             -- Run-level parameters and key statistics
#   02_uncertainty_by_feature.tsv  -- Per-feature metrics across all conditions
#   03_uncertainty_by_family.tsv   -- TE metrics aggregated per family
#   04_cv_summary_by_feature_type.tsv -- CV and InfRV summary: TEs vs genes
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. Install and load packages
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("dplyr", "matrixStats")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

bioc_packages <- c("tximport", "fishpond", "SummarizedExperiment")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
}

library(tximport)
library(fishpond)
library(SummarizedExperiment)
library(dplyr)
library(matrixStats)

# Prevent vroom/readr from spawning multiple threads.
# EAGAIN ("Resource temporarily unavailable") occurs when the system is near
# its thread or file-descriptor limit and a pthread_create() call fails.
# readr >= 2.x uses vroom internally; both the R option AND the environment
# variable must be set because they are independent controls.
options(readr.num_threads = 1)
Sys.setenv(VROOM_THREADS = "1")

# ------------------------------------------------------------------------------
# 1. User-defined parameters
# ------------------------------------------------------------------------------

# Path to the conditions TSV (columns: sample, condition)
CONDITIONS_FILE <- "~/rnaseq_d.tsv"

# Root directory containing one sub-folder per sample with quant.sf inside
SALMON_DIR <- "~/salmon_results"

# Directory for all output tables
OUT_DIR <- "~/te_uncertainty_results"

# TEdistill *.mapids classification file (no header, tab-separated)
CLASSIFICATION_FILE <- "~/Sviridis.flTE.mapids"

# Desired display order for conditions (must match labels in CONDITIONS_FILE)
# Adjust or extend this vector to fit your experimental design.
CONDITION_ORDER <- c("control", "drought", "rewatering")

# Pseudo-count added to Gibbs replicate values before metric calculation.
# Prevents division by zero (CV, phi) when a feature has zero counts.
# Common values: 1 (default) or 0.5.
PSEUDOCOUNT <- 1


# ------------------------------------------------------------------------------
# 2. Validate inputs
# ------------------------------------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("QUANTIFICATION UNCERTAINTY ANALYSIS\n")
cat(rep("=", 70), "\n\n", sep = "")

for (path in c(CONDITIONS_FILE, CLASSIFICATION_FILE)) {
  if (!file.exists(path)) stop("Required file not found: ", path)
}

if (!dir.exists(SALMON_DIR)) stop("Salmon directory not found: ", SALMON_DIR)

dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 3. Load conditions file
# ------------------------------------------------------------------------------

conditions <- read.table(
  CONDITIONS_FILE,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# Verify required columns
required_cols <- c("sample", "condition")
missing_cols  <- setdiff(required_cols, colnames(conditions))
if (length(missing_cols) > 0) {
  stop("CONDITIONS_FILE is missing column(s): ",
       paste(missing_cols, collapse = ", "))
}

cat("Conditions file loaded\n")
cat("  Samples:    ", nrow(conditions), "\n")
cat("  Conditions: ", paste(unique(conditions$condition), collapse = ", "), "\n\n")

# Normalise condition labels: trim whitespace then lowercase.
# Rows with empty condition values (after trimming) are dropped with a warning
# because they cannot be mapped to any condition and would silently propagate
# as "" keys throughout the pipeline, causing NULL lookups later.
conditions$condition <- tolower(trimws(conditions$condition))
CONDITION_ORDER      <- tolower(trimws(CONDITION_ORDER))

empty_rows <- nchar(conditions$condition) == 0 | is.na(conditions$condition)
if (any(empty_rows)) {
  warning(sum(empty_rows), " row(s) in CONDITIONS_FILE have an empty condition ",
          "value and will be ignored: ",
          paste(conditions$sample[empty_rows], collapse = ", "))
  conditions <- conditions[!empty_rows, ]
}

cat("  Conditions (normalised): ",
    paste(unique(conditions$condition), collapse = ", "), "\n\n")

# Derive condition levels from data, respecting CONDITION_ORDER where possible
detected_conditions <- unique(conditions$condition)
ordered_conditions  <- c(
  intersect(CONDITION_ORDER, detected_conditions),   # known, in order
  setdiff(detected_conditions, CONDITION_ORDER)      # any extras, appended
)


# ------------------------------------------------------------------------------
# 4. Locate Salmon quant.sf files
# ------------------------------------------------------------------------------

conditions$files <- file.path(SALMON_DIR, conditions$sample, "quant.sf")
missing_files    <- !file.exists(conditions$files)

if (any(missing_files)) {
  stop("quant.sf file(s) not found:\n",
       paste(conditions$files[missing_files], collapse = "\n"))
}

cat("All quant.sf files located\n\n")
files <- setNames(conditions$files, conditions$sample)


# ------------------------------------------------------------------------------
# 5. Probe first sample — get feature IDs only (no infReps loaded)
# ------------------------------------------------------------------------------
# Gibbs replicate matrices can exceed 1 GB per sample, making it impossible to
# load them via tximport even one sample at a time.  Instead we:
#   (a) Load counts/TPM only (dropInfReps=TRUE) to get feature IDs.
#   (b) Verify that the binary bootstraps.gz file exists for the first sample.
#   (c) Read Gibbs replicates directly in the per-sample loop (section 8),
#       one replicate at a time using readBin, accumulating running statistics.

cat("Probing first sample (feature IDs only, no infReps)...\n")
cat("  Reading: ", files[1], "\n", sep = "")
txi_first <- tryCatch(
  tximport(files[1], type = "salmon", txOut = TRUE, dropInfReps = TRUE),
  error = function(e) {
    cat("\n  ERROR: Failed to read first sample!\n")
    cat("    File:  ", files[1], "\n", sep = "")
    cat("    Reason: ", conditionMessage(e), "\n\n", sep = "")
    cat("  Checking file validity...\n")
    cat("    File exists:  ", file.exists(files[1]), "\n", sep = "")
    cat("    File size:    ", file.size(files[1]), " bytes\n", sep = "")
    cat("    First 500 chars:\n")
    try({
      con <- file(files[1], "r")
      lines <- readLines(con, n = 5)
      close(con)
      cat(paste("      ", lines), sep = "\n")
    }, silent = TRUE)
    stop("Cannot proceed without a valid first sample.")
  }
)

all_ids_first <- rownames(txi_first$counts)
n_all_first   <- length(all_ids_first)
rm(txi_first); gc()

# Locate the Salmon Gibbs binary file.  tximport stores both bootstrap and
# Gibbs samples at aux_info/bootstrap/bootstraps.gz regardless of samp_type.
boot_path_1 <- file.path(dirname(files[1]), "aux_info", "bootstrap", "bootstraps.gz")
if (!file.exists(boot_path_1)) {
  stop("Gibbs/bootstrap binary not found: ", boot_path_1,
       "\n  Re-run Salmon with --numGibbsSamples <N>.")
}

# Read replicate count from meta_info.json (tiny JSON, no memory impact)
minfo_1 <- jsonlite::fromJSON(
  file.path(dirname(files[1]), "aux_info", "meta_info.json")
)
n_gibbs <- minfo_1$num_bootstraps   # Salmon uses this field for both Gibbs and bootstrap
rm(minfo_1)

cat("  Features in first sample: ", n_all_first, "\n")
cat("  Gibbs replicates:         ", n_gibbs,      "\n")
cat("  bootstraps.gz confirmed:  ", boot_path_1, "\n\n")

# Metadata — derived from the conditions table (no SE needed)
samples           <- conditions$sample
n_samples         <- length(samples)
conditions_vector <- conditions$condition   # already tolower/trimws from section 3


# ------------------------------------------------------------------------------
# 6. Compute size factors (counts only, one sample at a time)
# ------------------------------------------------------------------------------
# Peak memory per iteration: one tximport object (~2 MB for a single quant.sf
# without infReps) + growing log_count_sum vector (n_features × 8 bytes).
# The count vectors for all samples are kept in a list (~n_features × n_samples
# × 8 bytes total) so each sample's size factor can be computed without a
# second file-read pass.

cat("Computing size factors (two-pass algorithm to minimize memory)...\n")
cat("  Pass 1: Accumulate log-counts to compute geometric means\n\n")

log_count_sum <- numeric(n_all_first)
samples_ok    <- logical(n_samples)

# PASS 1: Accumulate sums for geometric mean computation
for (s in seq_len(n_samples)) {
  qt <- tryCatch({
    read.table(files[s], header = TRUE, sep = "\t", nrows = n_all_first,
               colClasses = c("character", rep("numeric", 4)))
  }, error = function(e) {
    warning("Failed to read count data from ", files[s], ": ",
            conditionMessage(e))
    NULL
  })

  if (is.null(qt)) {
    samples_ok[s] <- FALSE
    next
  }

  if (!("NumReads" %in% colnames(qt))) {
    warning("NumReads column not found in ", files[s])
    samples_ok[s] <- FALSE
    next
  }

  samples_ok[s] <- TRUE
  cv            <- qt$NumReads
  rm(qt); gc()
  log_count_sum <- log_count_sum + log(cv + 0.5)
}

n_samples_ok <- sum(samples_ok)
cat("  Valid samples: ", n_samples_ok, " / ", n_samples, "\n", sep = "")

if (n_samples_ok < 2) {
  stop("Fewer than 2 valid samples; cannot compute size factors.")
}

geom_means_log <- log_count_sum / n_samples_ok
rm(log_count_sum); gc()

# PASS 2: Compute individual size factors without storing all counts
cat("  Pass 2: Compute size factors per sample\n\n")

size_factors <- numeric(n_samples)

for (s in seq_len(n_samples)) {
  if (!samples_ok[s]) {
    size_factors[s] <- NA_real_
    next
  }

  qt <- tryCatch({
    read.table(files[s], header = TRUE, sep = "\t", nrows = n_all_first,
               colClasses = c("character", rep("numeric", 4)))
  }, error = function(e) NULL)

  if (is.null(qt) || !("NumReads" %in% colnames(qt))) {
    size_factors[s] <- NA_real_
    next
  }

  cv <- qt$NumReads
  rm(qt); gc()
  size_factors[s] <- median(exp(log(cv + 0.5) - geom_means_log), na.rm = TRUE)
}
names(size_factors) <- samples
rm(geom_means_log); gc()

cat("  Size factors range: [",
    round(min(size_factors, na.rm = TRUE), 3), ", ",
    round(max(size_factors, na.rm = TRUE), 3), "]\n\n", sep = "")


# ------------------------------------------------------------------------------
# 7. Assign feature classification
# ------------------------------------------------------------------------------
# TE IDs follow the pattern: TE_<numericID>_copy<N>|Chr_XX:start-end|strand
# Gene IDs do NOT start with "TE_".
#
# Classification file format (TEdistill *.mapids, no header, tab-separated):
#   col 1: TE_<numericID>          (base ID, no copy suffix)
#   col 2: Family/Superfamily      (used here)
#   col 3: TE_<numericID>#Family/Superfamily
#
# All features are retained.  Genes receive feature_type = "gene" and NA for
# family/superfamily/classification.  TE copies are matched to the *.mapids
# file via their base ID (everything before "_copy"); unmatched copies receive
# classification = "Unknown".

te_class_ref <- read.table(
  CLASSIFICATION_FILE,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  col.names = c("te_id", "classification", "te_id_full")
)

cat("Classification file loaded: ", nrow(te_class_ref), " entries\n\n", sep = "")

all_ids   <- all_ids_first
n_all     <- length(all_ids)
is_te     <- grepl("^TE_", all_ids)
te_ids    <- all_ids[is_te]
gene_ids  <- all_ids[!is_te]
n_te      <- sum(is_te)
n_gene    <- sum(!is_te)

cat("Feature summary\n")
cat("  Total:  ", n_all,  "\n")
cat("  TEs:    ", n_te,   "\n")
cat("  Genes:  ", n_gene, "\n\n")

# TE classification lookup via base ID
te_base_id    <- sub("_copy.*", "", te_ids)
class_lookup  <- setNames(te_class_ref$classification, te_class_ref$te_id)
te_class_vals <- class_lookup[te_base_id]
te_class_vals[is.na(te_class_vals)]            <- "Unknown"
te_class_vals[te_class_vals == "unknown"]      <- "Unknown"

n_matched   <- sum(te_base_id %in% te_class_ref$te_id)
n_unmatched <- sum(!te_base_id %in% te_class_ref$te_id)
pct_matched <- round(n_matched / n_te * 100, 1)

cat("TE classification match rate\n")
cat("  Matched:   ", n_matched,   " (", pct_matched, "%)\n", sep = "")
cat("  Unmatched: ", n_unmatched, "\n\n")

te_family <- sub("/.*", "", te_class_vals)
te_family[te_family == "unknown"] <- "Unknown"

te_superfamily <- ifelse(
  grepl("/", te_class_vals),
  sub(".*/", "", te_class_vals),
  te_family
)
te_superfamily[te_superfamily == "unknown"] <- "Unknown"

# Build a classification data frame covering ALL features.
# Genes receive "gene" for classification/family/superfamily so they behave as
# a single coherent group in aggregated tables, distinguishable from
# unclassified TEs ("Unknown").
feature_class <- data.frame(
  feature_id     = all_ids,
  feature_type   = ifelse(is_te, "TE", "gene"),
  classification = c(te_class_vals, rep("gene", n_gene)),
  family         = c(te_family,     rep("gene", n_gene)),
  superfamily    = c(te_superfamily, rep("gene", n_gene)),
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 8. Compute per-sample uncertainty metrics (all features)
# ------------------------------------------------------------------------------
# For each feature × sample, five metrics summarise variation across the
# Gibbs replicates (count per sample determined during streaming, see n_gibbs).
#
# CV  (Coefficient of Variation) = SD / mean
#   Relative uncertainty; enables comparison across expression levels.
#
# SD  (Standard Deviation)
#   Absolute spread of replicates.
#
# IQR (Interquartile Range) = Q75 − Q25
#   Robust spread; less sensitive to extreme replicates.
#
# phi (inferential overdispersion) = variance / mean
#   Moment estimator of mapping ambiguity, as in Baldoni et al.
#   (Nucleic Acids Research, 2024, 52(3):e13).
#   Note: after scaleInfReps() the replicates are depth-adjusted, so phi
#   reflects scaled inferential overdispersion rather than quasi-Poisson
#   dispersion in the strict count-model sense.
#   rowVars() uses Bessel's correction (divides by n − 1).
#
# InfRV (Inferential Relative Variance) = max(var / mean² − 1/mean, 0)
#   Canonical fishpond metric introduced by Zhu et al.
#   (Nucleic Acids Research, 2019, 47(18):e105).  Subtracts the Poisson
#   sampling component (1/mean) from the squared CV, isolating the excess
#   variance due to read-assignment ambiguity.  More stable than CV for
#   features with moderate-to-high expression.  Computed via
#   fishpond::computeInfRV() on the full SummarizedExperiment and then
#   extracted per sample.
#
# Pseudo-count (PSEUDOCOUNT):
#   Added to all replicate values before metric calculation to avoid division
#   by zero in CV (SD / mean) and phi (variance / mean) when a feature has
#   zero counts.  Configured in section 1; default = 1.
#   Note: InfRV is computed by fishpond before the pseudo-count is applied,
#   using its own internal stabilisation.

cat("Computing uncertainty metrics for all features (one sample at a time)...\n")
cat("  Memory strategy: running condition sums instead of",
    "full n_features × n_samples matrices\n\n")

feat_names  <- all_ids
n_cond      <- length(ordered_conditions)

# Map each sample to its condition index for O(1) lookup in the loop
sample_cond <- match(conditions_vector, ordered_conditions)

# Allocate sum and non-NA-count matrices: n_all × n_cond.
# Peak memory ≈ n_all × n_cond × 8 bytes per matrix — independent of n_samples.
mk_sum <- function() matrix(0,  nrow = n_all, ncol = n_cond,
                             dimnames = list(feat_names, ordered_conditions))
mk_cnt <- function() matrix(0L, nrow = n_all, ncol = n_cond,
                             dimnames = list(feat_names, ordered_conditions))

tpm_s <- mk_sum(); tpm_n <- mk_cnt()
cv_s  <- mk_sum(); cv_n  <- mk_cnt()
sd_s  <- mk_sum(); sd_n  <- mk_cnt()
iqr_s <- mk_sum(); iqr_n <- mk_cnt()
phi_s <- mk_sum(); phi_n <- mk_cnt()
inf_s <- mk_sum(); inf_n <- mk_cnt()

# Global running sums across ALL samples (used for per-feature mean CV ranking
# and the run-level summary).  Vectors of length n_all.
cv_g_s  <- numeric(n_all)
cv_g_n  <- integer(n_all)
phi_g_s <- numeric(n_all)
inf_g_s <- numeric(n_all)

for (s in seq_len(n_samples)) {
  # Skip if this sample failed size factor computation
  if (!samples_ok[s]) {
    cat("  Sample ", s, "/", n_samples, ": ", samples[s],
        " [SKIPPED — invalid size factors]\n", sep = "")
    next
  }

  ci <- sample_cond[s]
  cat("  Sample ", s, "/", n_samples,
      " [", ordered_conditions[ci], "]: ", samples[s], "\n", sep = "")

  # ── TPM (reads directly from quant.sf, avoiding tximport overhead) ──────────
  qt <- tryCatch({
    read.table(files[s], header = TRUE, sep = "\t", nrows = n_all_first,
               colClasses = c("character", rep("numeric", 4)))
  }, error = function(e) {
    warning("Failed to read TPM from ", files[s], ": ", conditionMessage(e))
    NULL
  })

  if (is.null(qt) || !("TPM" %in% colnames(qt))) {
    warning("  Sample ", samples[s], " will be skipped entirely (missing TPM).")
    next
  }

  tpm_v <- qt$TPM
  rm(qt); gc()

  ok    <- !is.na(tpm_v)
  tpm_s[ok, ci] <- tpm_s[ok, ci] + tpm_v[ok]
  tpm_n[ok, ci] <- tpm_n[ok, ci] + 1L
  rm(tpm_v)

  # ── Gibbs replicates: stream one replicate at a time from binary file ──────
  # Salmon stores all replicates in aux_info/bootstrap/bootstraps.gz as a flat
  # sequence of float64 values written column-major: first n_features values =
  # replicate 1, next n_features values = replicate 2, etc.
  # Reading n_features doubles repeatedly with readBin streams exactly one
  # replicate per call, keeping peak memory at ~4 × n_features × 8 bytes.
  boot_path <- file.path(dirname(files[s]), "aux_info", "bootstrap",
                         "bootstraps.gz")

  if (!file.exists(boot_path)) {
    warning("No bootstraps.gz for sample ", samples[s], "; skipping metrics.")
    next
  }

  sf_inv   <- 1 / size_factors[s]
  sum_raw  <- numeric(n_all)   # scaled values — used for InfRV
  sum_raw2 <- numeric(n_all)
  sum_pc   <- numeric(n_all)   # scaled + pseudocount — used for CV/SD/phi
  sum_pc2  <- numeric(n_all)
  n_reps   <- 0L

  boot_con <- NULL
  tryCatch({
    boot_con <- gzcon(file(boot_path, "rb"))
    repeat {
      vals <- readBin(boot_con, what = "double", n = n_all, size = 8L,
                      endian = "little")
      if (length(vals) == 0L) break
      if (length(vals) != n_all) {
        warning("Partial replicate at end of ", boot_path, "; discarded.")
        break
      }
      vals_sc  <- vals * sf_inv
      vals_pc  <- vals_sc + PSEUDOCOUNT
      sum_raw  <- sum_raw  + vals_sc
      sum_raw2 <- sum_raw2 + vals_sc  * vals_sc
      sum_pc   <- sum_pc   + vals_pc
      sum_pc2  <- sum_pc2  + vals_pc  * vals_pc
      n_reps   <- n_reps   + 1L
    }
  }, error = function(e) {
    warning("Error reading bootstraps.gz for sample ", samples[s],
            "\n  Error: ", conditionMessage(e))
  }, finally = {
    if (!is.null(boot_con)) close(boot_con)
  })

  if (n_reps < 2L) {
    warning("Fewer than 2 replicates for sample ", samples[s],
            "; uncertainty metrics skipped.")
    next
  }

  # ── Compute metrics from running sums (Bessel-corrected variance) ──────────
  mu_raw  <- sum_raw / n_reps
  mu_pc   <- sum_pc  / n_reps
  var_raw <- pmax((sum_raw2 - n_reps * mu_raw^2) / (n_reps - 1L), 0)
  var_pc  <- pmax((sum_pc2  - n_reps * mu_pc^2)  / (n_reps - 1L), 0)
  rm(sum_raw, sum_raw2, sum_pc, sum_pc2)

  sd_v  <- sqrt(var_pc)
  cv_v  <- sd_v  / mu_pc
  iqr_v <- 1.35 * sd_v   # normal-distribution approximation: IQR ≈ 1.35 SD
  phi_v <- var_pc / mu_pc

  infrv <- pmax(var_raw / mu_raw^2 - 1 / mu_raw, 0)
  infrv[!is.finite(infrv)] <- NA
  rm(var_raw, var_pc, mu_raw, mu_pc); gc()

  ok_i  <- !is.na(infrv)
  inf_s[ok_i, ci] <- inf_s[ok_i, ci] + infrv[ok_i]
  inf_n[ok_i, ci] <- inf_n[ok_i, ci] + 1L
  inf_g_s         <- inf_g_s + ifelse(ok_i, infrv, 0)

  ok_cv <- !is.na(cv_v)
  cv_s[ok_cv, ci] <- cv_s[ok_cv, ci] + cv_v[ok_cv]
  cv_n[ok_cv, ci] <- cv_n[ok_cv, ci] + 1L
  cv_g_s[ok_cv]   <- cv_g_s[ok_cv] + cv_v[ok_cv]
  cv_g_n[ok_cv]   <- cv_g_n[ok_cv] + 1L

  ok_sd <- !is.na(sd_v)
  sd_s[ok_sd, ci] <- sd_s[ok_sd, ci] + sd_v[ok_sd]
  sd_n[ok_sd, ci] <- sd_n[ok_sd, ci] + 1L

  ok_iq <- !is.na(iqr_v)
  iqr_s[ok_iq, ci] <- iqr_s[ok_iq, ci] + iqr_v[ok_iq]
  iqr_n[ok_iq, ci] <- iqr_n[ok_iq, ci] + 1L

  ok_ph <- !is.na(phi_v)
  phi_s[ok_ph, ci] <- phi_s[ok_ph, ci] + phi_v[ok_ph]
  phi_n[ok_ph, ci] <- phi_n[ok_ph, ci] + 1L
  phi_g_s          <- phi_g_s + ifelse(ok_ph, phi_v, 0)
}

gc()
cat("\nMetrics computed\n\n")


# ------------------------------------------------------------------------------
# 9. Sanitise feature ID names
# ------------------------------------------------------------------------------
# Some TE IDs contain "|" (e.g. "TE_001|Chr01:100-200|+").
# Retain only the portion before the first "|".

clean_feat_names <- sub("\\|.*", "", feat_names)

rownames(tpm_s) <- rownames(tpm_n) <- clean_feat_names
rownames(cv_s)  <- rownames(cv_n)  <- clean_feat_names
rownames(sd_s)  <- rownames(sd_n)  <- clean_feat_names
rownames(iqr_s) <- rownames(iqr_n) <- clean_feat_names
rownames(phi_s) <- rownames(phi_n) <- clean_feat_names
rownames(inf_s) <- rownames(inf_n) <- clean_feat_names
names(cv_g_s)   <- names(cv_g_n)   <- clean_feat_names
names(phi_g_s)  <- names(inf_g_s)  <- clean_feat_names
feat_names      <- clean_feat_names
feature_class$feature_id <- sub("\\|.*", "", feature_class$feature_id)

cat("Feature IDs sanitised (\"|\"-truncated)\n")
cat("  TE example:   ", head(feat_names[is_te],  1), "\n", sep = "")
cat("  Gene example: ", head(feat_names[!is_te], 1), "\n\n", sep = "")

# Global mean CV per feature across all samples (used for ranking)
global_cv_all  <- cv_g_s / ifelse(cv_g_n > 0L, cv_g_n, NA_integer_)
global_cv_te   <- global_cv_all[is_te]
global_cv_gene <- global_cv_all[!is_te]


# ------------------------------------------------------------------------------
# 10. Compute condition means and build _by_cond lists
# ------------------------------------------------------------------------------
# Divide accumulated sums by counts; where count == 0 the result is NA.
# Free each pair of sum/count matrices immediately after use.

safe_div <- function(s_mat, n_mat) {
  res           <- s_mat / n_mat
  res[n_mat == 0L] <- NA_real_
  res
}

# Helper: turn an n_all × n_cond mean matrix into the named list expected by
# build_condition_columns.  unname() removes rowname artifacts from column
# extraction so no "|" characters leak into the output data frame.
mat_to_list <- function(mat) {
  lst <- lapply(ordered_conditions, function(cond) unname(mat[, cond]))
  names(lst) <- ordered_conditions
  lst
}

tpm_mean      <- safe_div(tpm_s, tpm_n); rm(tpm_s, tpm_n)
tpm_by_cond   <- mat_to_list(tpm_mean);  rm(tpm_mean)

cv_mean       <- safe_div(cv_s,  cv_n);  rm(cv_s,  cv_n)
cv_by_cond    <- mat_to_list(cv_mean);   rm(cv_mean)

sd_mean       <- safe_div(sd_s,  sd_n);  rm(sd_s,  sd_n)
sd_by_cond    <- mat_to_list(sd_mean);   rm(sd_mean)

iqr_mean      <- safe_div(iqr_s, iqr_n); rm(iqr_s, iqr_n)
iqr_by_cond   <- mat_to_list(iqr_mean);  rm(iqr_mean)

phi_mean      <- safe_div(phi_s, phi_n); rm(phi_s, phi_n)
phi_by_cond   <- mat_to_list(phi_mean);  rm(phi_mean)

infrv_mean    <- safe_div(inf_s, inf_n); rm(inf_s, inf_n)
infrv_by_cond <- mat_to_list(infrv_mean); rm(infrv_mean)

gc()


# ------------------------------------------------------------------------------
# 11. Build per-feature metrics table
# ------------------------------------------------------------------------------

# Dynamically create one column set per condition.
# metric_by_cond must be a named list (names = ordered_conditions).
# If a condition key is missing, a column of NA is inserted and a warning
# is emitted so the problem is diagnosable without crashing the script.
build_condition_columns <- function(metric_by_cond, prefix) {
  cols <- lapply(ordered_conditions, function(cond) {
    v <- metric_by_cond[[cond]]
    if (is.null(v)) {
      warning("build_condition_columns: '", prefix, "[[", cond, "]]' is NULL ",
              "(available names: ",
              paste(names(metric_by_cond), collapse = ", "), "). ",
              "Filling with NA.")
      return(rep(NA_real_, n_all))
    }
    if (length(v) != n_all) {
      warning("build_condition_columns: '", prefix, "[[", cond, "]]' has ",
              "length ", length(v), " (expected ", n_all, ").")
    }
    v
  })
  names(cols) <- paste0(prefix, "_", ordered_conditions)
  cols
}

copy_metrics <- as.data.frame(
  c(
    list(feature_id = feat_names),
    build_condition_columns(tpm_by_cond,   "mean_tpm"),
    build_condition_columns(cv_by_cond,    "mean_cv"),
    build_condition_columns(sd_by_cond,    "mean_sd"),
    build_condition_columns(iqr_by_cond,   "mean_iqr"),
    build_condition_columns(phi_by_cond,   "mean_phi"),
    build_condition_columns(infrv_by_cond, "mean_infrv")
  ),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
rownames(copy_metrics) <- NULL

# Global metrics across all conditions
cv_cols    <- paste0("mean_cv_",    ordered_conditions)
phi_cols   <- paste0("mean_phi_",   ordered_conditions)
infrv_cols <- paste0("mean_infrv_", ordered_conditions)

copy_metrics$mean_cv_global <- rowMeans(copy_metrics[, cv_cols, drop = FALSE],
                                        na.rm = TRUE)
copy_metrics$mean_phi_global <- rowMeans(copy_metrics[, phi_cols, drop = FALSE],
                                         na.rm = TRUE)
copy_metrics$mean_infrv_global <- rowMeans(copy_metrics[, infrv_cols, drop = FALSE],
                                           na.rm = TRUE)
copy_metrics$mean_tpm_global <- rowMeans(
  copy_metrics[, paste0("mean_tpm_", ordered_conditions), drop = FALSE],
  na.rm = TRUE
)

# Attach classification (feature_type, family, superfamily, classification)
copy_metrics <- merge(
  copy_metrics,
  feature_class[, c("feature_id", "feature_type",
                    "family", "superfamily", "classification")],
  by    = "feature_id",
  all.x = TRUE
)

# Sort by global CV (descending) and add rank
copy_metrics <- copy_metrics[
  order(copy_metrics$mean_cv_global, decreasing = TRUE, na.last = TRUE), ]
copy_metrics$cv_rank <- seq_len(nrow(copy_metrics))

write.table(
  copy_metrics,
  file.path(OUT_DIR, "tables", "02_uncertainty_by_feature.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Table saved: 02_uncertainty_by_feature.tsv\n")
cat("  Most uncertain feature: ", copy_metrics$feature_id[1],
    "  (global CV = ", round(copy_metrics$mean_cv_global[1], 4), ")\n\n",
    sep = "")


# ------------------------------------------------------------------------------
# 14. Aggregate metrics by TE family
# ------------------------------------------------------------------------------

family_metrics <- copy_metrics %>%
  group_by(family, superfamily) %>%
  summarise(
    n_copies = n(),
    mean_cv_global = mean(mean_cv_global, na.rm = TRUE),
    median_cv_global = median(mean_cv_global, na.rm = TRUE),
    mean_phi_global = mean(mean_phi_global, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_cv_global))

write.table(
  family_metrics,
  file.path(OUT_DIR, "tables", "03_uncertainty_by_family.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Table saved: 03_uncertainty_by_family.tsv\n\n")


# ------------------------------------------------------------------------------
# 13. CV and InfRV summary by feature type (TE vs gene)
# ------------------------------------------------------------------------------
# Summarises uncertainty metrics separately for TEs and genes, providing a
# direct quantitative comparison of quantification reliability between the
# two feature types across all samples and conditions.

cv_summary <- do.call(rbind, lapply(c("TE", "gene"), function(ftype) {
  idx  <- copy_metrics$feature_type == ftype
  vals_cv    <- copy_metrics$mean_cv_global[idx]
  vals_infrv <- copy_metrics$mean_infrv_global[idx]
  data.frame(
    feature_type    = ftype,
    n_features      = sum(!is.na(vals_cv)),
    mean_cv         = round(mean(vals_cv,    na.rm = TRUE), 4),
    median_cv       = round(median(vals_cv,  na.rm = TRUE), 4),
    sd_cv           = round(sd(vals_cv,      na.rm = TRUE), 4),
    p25_cv          = round(quantile(vals_cv,    0.25, na.rm = TRUE), 4),
    p75_cv          = round(quantile(vals_cv,    0.75, na.rm = TRUE), 4),
    mean_infrv      = round(mean(vals_infrv,   na.rm = TRUE), 4),
    median_infrv    = round(median(vals_infrv, na.rm = TRUE), 4),
    stringsAsFactors = FALSE
  )
}))

write.table(
  cv_summary,
  file.path(OUT_DIR, "tables", "04_cv_summary_by_feature_type.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Table saved: 04_cv_summary_by_feature_type.tsv\n\n")


# ------------------------------------------------------------------------------
# 16. Run-level summary table
# ------------------------------------------------------------------------------

n_complete_te <- sum(cv_g_n == n_samples)
n_families    <- length(unique(copy_metrics$family[
  copy_metrics$family != "Unknown"]))

run_summary <- data.frame(
  parameter = c(
    "conditions_file", "salmon_dir", "classification_file",
    "n_gibbs_replicates_sample1", "n_samples", "n_conditions",
    "n_te_copies", "n_genes",
    "n_te_families", "te_classification_match_pct",
    "n_te_complete_cv", "pct_te_complete_cv",
    "global_mean_cv_te", "global_mean_cv_genes",
    "most_uncertain_feature", "most_uncertain_feature_cv"
  ),
  value = c(
    CONDITIONS_FILE, SALMON_DIR, CLASSIFICATION_FILE,
    n_gibbs, n_samples, length(ordered_conditions),
    n_te, n_gene,
    n_families, pct_matched,
    n_complete_te, round(n_complete_te / n_te * 100, 1),
    round(mean(global_cv_te,   na.rm = TRUE), 4),
    round(mean(global_cv_gene, na.rm = TRUE), 4),
    copy_metrics$feature_id[1],
    round(copy_metrics$mean_cv_global[1], 4)
  ),
  stringsAsFactors = FALSE
)

write.table(
  run_summary,
  file.path(OUT_DIR, "tables", "01_run_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("Table saved: 01_run_summary.tsv\n\n")


# ------------------------------------------------------------------------------
# 17. Final summary to console
# ------------------------------------------------------------------------------

cat(rep("=", 70), "\n", sep = "")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nRun parameters\n")
cat("  Conditions file:    ", CONDITIONS_FILE,     "\n")
cat("  Salmon directory:   ", SALMON_DIR,          "\n")
cat("  Classification file:", CLASSIFICATION_FILE, "\n")
cat("  Output directory:   ", OUT_DIR,             "\n")

cat("\nDataset\n")
cat("  TE copies:           ", n_te,      "\n")
cat("  Protein-coding genes:", n_gene,    "\n")
cat("  Samples:             ", n_samples, "\n")
cat("  Gibbs replicates:    ", n_gibbs, " (from first sample meta_info.json)\n")
cat("  TE families:         ", n_families,"\n")
cat("  Classification match:", pct_matched, "%\n")
cat("  TEs with complete CV:",
    n_complete_te, " (", round(n_complete_te / n_te * 100, 1), "%)\n",
    sep = "")

cat("\nMean CV by condition\n")
for (cond in ordered_conditions) {
  col <- paste0("mean_cv_", cond)
  cat("  ", formatC(cond, width = 14, flag = "-"), ": ",
      round(mean(copy_metrics[[col]], na.rm = TRUE), 4), "\n", sep = "")
}

cat("\nGlobal mean CV\n")
cat("  TEs:   ", round(mean(global_cv_te,   na.rm = TRUE), 4), "\n")
cat("  Genes: ", round(mean(global_cv_gene, na.rm = TRUE), 4), "\n")

cat("\nMost uncertain feature\n")
cat("  ", copy_metrics$feature_id[1],
    "  (CV = ", round(copy_metrics$mean_cv_global[1], 4), ")\n\n", sep = "")

cat("Output tables (", file.path(OUT_DIR, "tables"), "/)\n", sep = "")
cat("  01_run_summary.tsv\n")
cat("  02_uncertainty_by_feature.tsv\n")
cat("  03_uncertainty_by_family.tsv\n")
cat("  04_cv_summary_by_feature_type.tsv\n")
cat(rep("=", 70), "\n\n", sep = "")

sessionInfo()
