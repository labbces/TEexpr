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

# Normalise condition labels to lowercase to avoid silent mismatches caused
# by capitalisation differences (e.g. "Control" vs "control")
conditions$condition <- tolower(conditions$condition)
CONDITION_ORDER      <- tolower(CONDITION_ORDER)

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
# 5. Import Salmon output with tximport
# ------------------------------------------------------------------------------

cat("Importing Salmon data (this may take several minutes)...\n")

txi <- tximport(
  files,
  type        = "salmon",
  txOut       = TRUE,
  dropInfReps = FALSE
)

n_gibbs <- ncol(txi$infReps[[1]])

cat("  Features imported:     ", nrow(txi$counts), "\n")
cat("  Samples detected:      ", ncol(txi$counts), "\n")
cat("  Gibbs replicates / feature: ", n_gibbs, "\n\n")

if (n_gibbs == 0) {
  stop("No Gibbs replicates found. ",
       "Re-run Salmon with --numGibbsSamples <N>.")
}


# ------------------------------------------------------------------------------
# 6. Build SummarizedExperiment
# ------------------------------------------------------------------------------
# tximport returns txi$infReps as a list of matrices (one per sample), each of
# dimension n_features × n_Gibbs.  We transpose the indexing so that each assay
# 'infRepK' is an n_features × n_samples matrix, consistent with tximport
# convention and required by fishpond::scaleInfReps().
#
# Resulting assays:
#   counts     -- read counts             [n_features × n_samples]
#   abundance  -- TPM values              [n_features × n_samples]
#   length     -- effective lengths       [n_features × n_samples]
#   infRep1 .. infRep<n_gibbs>            [n_features × n_samples]

cat("Restructuring inferential replicates for SummarizedExperiment...\n")

inf_assays <- setNames(
  lapply(seq_len(n_gibbs), function(k) {
    mat <- do.call(cbind, lapply(txi$infReps, function(s_mat) {
      s_mat[, k, drop = FALSE]
    }))
    colnames(mat) <- conditions$sample
    mat
  }),
  paste0("infRep", seq_len(n_gibbs))
)

se <- SummarizedExperiment(
  assays  = c(
    list(
      counts    = txi$counts,
      abundance = txi$abundance,
      length    = txi$length
    ),
    inf_assays
  ),
  colData = DataFrame(
    condition = factor(conditions$condition,
                       levels = ordered_conditions),
    row.names = conditions$sample
  )
)

# Scale inferential replicates to correct for differences in sequencing depth
se <- scaleInfReps(se)

samples           <- colnames(se)
n_samples         <- length(samples)
conditions_vector <- as.character(colData(se)$condition)

cat("  SummarizedExperiment built and replicates scaled\n")
cat("  Dimensions: ", nrow(se), " features × ", n_samples, " samples\n\n",
    sep = "")


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

all_ids   <- rownames(se)
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
# n_gibbs Gibbs replicates.
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

cat("Computing uncertainty metrics for all features...\n")

# --- InfRV via fishpond -------------------------------------------------------
# computeInfRV() adds "mean" and "variance" assays to se (averaged across
# Gibbs replicates per feature per sample) and stores InfRV in rowData(se).
# We compute the InfRV matrix directly from those two assays using the
# canonical formula: InfRV = max(variance / mean² - 1/mean, 0)
# This is equivalent to what fishpond computes internally and avoids any
# version-dependent naming of the output assay.
se <- computeInfRV(se)

infrv_mean <- assay(se, "mean")
infrv_var  <- assay(se, "variance")
mat_infrv  <- pmax(infrv_var / infrv_mean^2 - 1 / infrv_mean, 0)
mat_infrv[!is.finite(mat_infrv)] <- NA

rep_names  <- paste0("infRep", seq_len(n_gibbs))
feat_names <- all_ids

mat_cv   <- matrix(NA, nrow = n_all, ncol = n_samples,
                   dimnames = list(feat_names, samples))
mat_sd   <- mat_cv
mat_iqr  <- mat_cv
mat_phi  <- mat_cv
mat_mean <- mat_cv

for (s in seq_len(n_samples)) {
  cat("  Sample ", s, "/", n_samples, ": ", samples[s], "\n", sep = "")
  rep_mat   <- sapply(rep_names, function(r) assay(se, r)[, s]) + PSEUDOCOUNT
  mean_vals <- rowMeans(rep_mat, na.rm = TRUE)
  var_vals  <- matrixStats::rowVars(rep_mat, na.rm = TRUE)
  sd_vals   <- sqrt(var_vals)
  iqr_vals  <- matrixStats::rowIQRs(rep_mat, na.rm = TRUE)
  
  mat_mean[, s] <- mean_vals
  mat_cv[, s]   <- sd_vals / mean_vals
  mat_sd[, s]   <- sd_vals
  mat_iqr[, s]  <- iqr_vals
  mat_phi[, s]  <- var_vals / mean_vals
}

cat("\nMetrics computed\n\n")


# ------------------------------------------------------------------------------
# 9. Sanitise feature ID names
# ------------------------------------------------------------------------------
# Some TE IDs contain "|" (e.g. "TE_001|Chr01:100-200|+"), which can disrupt
# data.frame construction.  We retain only the portion before the first "|".
# Gene IDs typically do not contain "|", so this step is safe for all features.

clean_feat_names <- sub("\\|.*", "", rownames(mat_cv))

rownames(mat_cv)    <- clean_feat_names
rownames(mat_sd)    <- clean_feat_names
rownames(mat_iqr)   <- clean_feat_names
rownames(mat_phi)   <- clean_feat_names
rownames(mat_mean)  <- clean_feat_names
rownames(mat_infrv) <- clean_feat_names
feat_names          <- clean_feat_names
feature_class$feature_id <- sub("\\|.*", "", feature_class$feature_id)

cat("Feature IDs sanitised (\"|\"-truncated)\n")
cat("  TE example:   ", head(feat_names[is_te],  1), "\n", sep = "")
cat("  Gene example: ", head(feat_names[!is_te], 1), "\n\n", sep = "")

# Global mean CV per feature across all samples (used for ranking)
global_cv_all  <- rowMeans(mat_cv, na.rm = TRUE)
global_cv_te   <- global_cv_all[is_te]
global_cv_gene <- global_cv_all[!is_te]


# ------------------------------------------------------------------------------
# 10. Summarise metrics by condition
# ------------------------------------------------------------------------------

# Helper: column-wise mean for samples belonging to a given condition
summarise_by_condition <- function(mat, cond) {
  cols <- which(conditions_vector == cond)
  if (length(cols) == 1) return(mat[, cols])
  rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
}

cv_by_cond <- lapply(ordered_conditions,
                     function(c) summarise_by_condition(mat_cv, c))
sd_by_cond <- lapply(ordered_conditions,
                     function(c) summarise_by_condition(mat_sd, c))
iqr_by_cond <- lapply(ordered_conditions,
                      function(c) summarise_by_condition(mat_iqr, c))
phi_by_cond <- lapply(ordered_conditions,
                      function(c) summarise_by_condition(mat_phi, c))
infrv_by_cond <- lapply(ordered_conditions,
                        function(c) summarise_by_condition(mat_infrv, c))

tpm_mat     <- assay(se, "abundance")
tpm_by_cond <- lapply(ordered_conditions, function(c) {
  cols <- which(conditions_vector == c)
  if (length(cols) == 1) return(tpm_mat[, cols])
  rowMeans(tpm_mat[, cols, drop = FALSE], na.rm = TRUE)
})

names(cv_by_cond) <- names(sd_by_cond)    <-
  names(iqr_by_cond) <- names(phi_by_cond) <-
  names(infrv_by_cond) <- names(tpm_by_cond) <- ordered_conditions

# Drop vector names to prevent "|" characters from leaking into data frames
strip_names <- function(lst) lapply(lst, function(x) { names(x) <- NULL; x })
cv_by_cond    <- strip_names(cv_by_cond)
sd_by_cond    <- strip_names(sd_by_cond)
iqr_by_cond   <- strip_names(iqr_by_cond)
phi_by_cond   <- strip_names(phi_by_cond)
infrv_by_cond <- strip_names(infrv_by_cond)
tpm_by_cond   <- strip_names(tpm_by_cond)


# ------------------------------------------------------------------------------
# 11. Build per-feature metrics table
# ------------------------------------------------------------------------------

# Dynamically create one column set per condition
build_condition_columns <- function(metric_by_cond, prefix) {
  cols <- lapply(ordered_conditions, function(c) metric_by_cond[[c]])
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
  stringsAsFactors = FALSE
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

n_complete_te <- sum(complete.cases(mat_cv))
n_families    <- length(unique(copy_metrics$family[
  copy_metrics$family != "Unknown"]))

run_summary <- data.frame(
  parameter = c(
    "conditions_file", "salmon_dir", "classification_file",
    "n_gibbs_replicates", "n_samples", "n_conditions",
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
cat("  Gibbs replicates:    ", n_gibbs,   "\n")
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
