# ==============================================================================
# TE quantification uncertainty analysis using Salmon Gibbs replicates.
#
# Compatible with single and multiple samples.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Install and load packages
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("ggplot2", "dplyr", "tidyr", "matrixStats", "viridis")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_packages <- c("tximport", "fishpond", "SummarizedExperiment")

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(tximport)
library(fishpond)
library(SummarizedExperiment)
library(ggplot2)
library(dplyr)
library(tidyr)
library(matrixStats)
library(viridis)


# ------------------------------------------------------------------------------
# Set paths
# ------------------------------------------------------------------------------

CONDITIONS_FILE     <- "~/rnaseq_d.tsv"
SALMON_DIR          <- "~/salmon_results"
OUT_DIR             <- "~/test_results_newplots"
CLASSIFICATION_FILE <- "~/Sviridis.flTE.mapids"

# Number of top uncertain TEs to show in expression profile plot
N_TOP_TE_PROFILES <- 12

if (!file.exists(CONDITIONS_FILE)) {
  stop("Conditions file not found: ", CONDITIONS_FILE)
}

if (!file.exists(CLASSIFICATION_FILE)) {
  stop("Classification file not found: ", CLASSIFICATION_FILE)
}


# ------------------------------------------------------------------------------
# Read the conditions file
# ------------------------------------------------------------------------------

conditions <- read.table(
  CONDITIONS_FILE,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

cat("Conditions file loaded!\n")
cat("Number of samples:", nrow(conditions), "\n")
cat("Groups found:", paste(unique(conditions$condition), collapse = ", "), "\n\n")
dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------------------------
# Locate Salmon output files
# ------------------------------------------------------------------------------

conditions$files <- file.path(SALMON_DIR, conditions$sample, "quant.sf")

missing_files <- !file.exists(conditions$files)
if (any(missing_files)) {
  stop("The following quant.sf files were not found:\n",
       paste(conditions$files[missing_files], collapse = "\n"))
}

cat("All quant.sf files found!\n\n")

files <- setNames(conditions$files, conditions$sample)


# ------------------------------------------------------------------------------
# Import Salmon data with tximport
# ------------------------------------------------------------------------------

cat("Importing Salmon data (this may take a few minutes)...\n")

txi <- tximport(
  files,
  type        = "salmon",
  txOut       = TRUE,
  dropInfReps = FALSE
)

n_gibbs <- ncol(txi$infReps[[1]])

cat("TEs imported:         ", nrow(txi$counts), "\n")
cat("Samples detected:     ", ncol(txi$counts), "\n")
cat("Gibbs replicates/TE:  ", n_gibbs, "\n\n")

if (n_gibbs == 0) {
  stop("No Gibbs replicates found. Make sure Salmon was run with --numGibbsSamples.")
}


# ------------------------------------------------------------------------------
# Build a SummarizedExperiment object
# ------------------------------------------------------------------------------
# tximport returns txi$infReps as a list of 200 matrices, each with dimensions
# n_TEs x n_samples. Then they are assigned sequential names and stored as
# individual assays alongside counts, abundance, and length.
#
# The resulting object (se) centralizes all data with guaranteed row/column
# alignment across assays:
#   assay(se, "counts")    -> read counts       [n_TEs x n_samples]
#   assay(se, "abundance") -> TPM values        [n_TEs x n_samples]
#   assay(se, "infRep1")   -> Gibbs replicate 1 [n_TEs x n_samples]
#   assay(se, "infRep200") -> Gibbs replicate 200 [n_TEs x n_samples]
#   colData(se)$condition  -> condition label per sample

cat("Restructuring infReps for SummarizedExperiment...\n")

inf_assays <- setNames(
  lapply(seq_len(n_gibbs), function(k) {
    mat <- do.call(cbind, lapply(txi$infReps, function(sample_mat) {
      sample_mat[, k, drop = FALSE]
    }))
    colnames(mat) <- conditions$sample
    mat
  }),
  paste0("infRep", seq_len(n_gibbs))
)

cat("infRep1 after restructuring:", dim(inf_assays[[1]]),
    "(expected:", nrow(txi$counts), "x", ncol(txi$counts), ")\n\n")

se <- SummarizedExperiment(
  assays  = c(
    list(counts    = txi$counts,
         abundance = txi$abundance,
         length    = txi$length),
    inf_assays
  ),
  colData = DataFrame(
    condition = factor(conditions$condition),
    row.names = conditions$sample
  )
)

# Scale inferential replicates (The amount of reads in each of the
# conditions is different. It's required to scale this or do a correction for
# the sequencing effort).
se <- scaleInfReps(se)
samples   <- colnames(se)
n_samples <- length(samples)
# se now contains scaled inferential replicates

cat("SummarizedExperiment built!\n")


# ------------------------------------------------------------------------------
# Assign TE classification from external file
# ------------------------------------------------------------------------------
# TE IDs follow: TE_<numericID>_copy<N>|Chr_XX:start-end|strand
# Genes are IDs NOT starting with "TE_".
# TEdistill "[specie].mapids" file: no header, tab-separated, 3 columns:
#   col 1: TE_<numericID>
#   col 2: ClassificationA/a  <- used here
#   col 3: TE_<numericID>#ClassificationA/a

te_class_ref <- read.table(
  CLASSIFICATION_FILE,
  header           = FALSE,
  sep              = "\t",
  stringsAsFactors = FALSE,
  col.names        = c("te_id", "classification", "te_id_full")
)

cat("Classification file loaded:", nrow(te_class_ref), "entries\n\n")

all_ids  <- rownames(se)
is_te    <- grepl("^TE_", all_ids)
te_ids   <- all_ids[is_te]
gene_ids <- all_ids[!is_te]

cat("Total IDs:  ", length(all_ids),  "\n")
cat("TEs found:  ", length(te_ids),   "\n")
cat("Genes found:", length(gene_ids), "\n\n")

# Extract TE_<numericID> by removing everything from "_copy" onward
te_base_id <- sub("_copy.*", "", te_ids)

# Check coverage
n_matched   <- sum(te_base_id %in% te_class_ref$te_id)
n_unmatched <- sum(!te_base_id %in% te_class_ref$te_id)
cat("Matched:  ", n_matched,   "\n")
cat("Unmatched:", n_unmatched, "\n\n")

# Map classification using lookup
class_lookup      <- setNames(te_class_ref$classification, te_class_ref$te_id)
te_classification <- class_lookup[te_base_id]

# Fill unmatched
te_classification[is.na(te_classification)] <- "Unknown"

# Split into family and superfamily
te_family      <- sub("/.*", "", te_classification)
te_superfamily <- ifelse(
  grepl("/", te_classification),
  sub(".*/", "", te_classification),
  te_family
)

# Standardize ambiguous labels
te_classification[te_classification == "unknown"] <- "Unknown"
te_family[te_family                 == "unknown"] <- "Unknown"
te_superfamily[te_superfamily       == "unknown"] <- "Unknown"


# ------------------------------------------------------------------------------
# Compute gene CV before filtering SE (needed for Plot 01)
# ------------------------------------------------------------------------------
# Plot 01 compares CV distributions of TEs vs protein-coding genes.
# Gene rows are only accessible before se is filtered to TEs, so we compute
# their global mean CV here and store it as global_cv_gene.

cat("Computing CV for protein-coding genes (for Plot 01)...\n")

rep_names <- paste0("infRep", seq_len(n_gibbs))
gene_idx  <- which(!is_te)
n_gene    <- length(gene_idx)

mat_cv_gene <- matrix(NA, nrow = n_gene, ncol = n_samples,
                      dimnames = list(gene_ids, samples))

for (s in seq_len(n_samples)) {
  rep_mat_g        <- sapply(rep_names, function(r) assay(se, r)[gene_idx, s]) + 1
  mean_g           <- rowMeans(rep_mat_g, na.rm = TRUE)
  sd_g             <- matrixStats::rowSds(rep_mat_g, na.rm = TRUE)
  mat_cv_gene[, s] <- sd_g / mean_g
}

global_cv_gene <- rowMeans(mat_cv_gene, na.rm = TRUE)
cat("  Done:", n_gene, "genes processed\n\n")


# ------------------------------------------------------------------------------
# Filter SE to retain only TE rows
# ------------------------------------------------------------------------------

se       <- se[is_te, ]
te_names <- rownames(se)
n_te     <- nrow(se)

cat("SE filtered to TEs only.\n")
cat("Dimensions:", n_te, "TEs x", n_samples, "samples\n\n")

if (length(te_family) != nrow(se) || length(te_superfamily) != nrow(se)) {
  stop("Alignment error: te_family or te_superfamily length does not match nrow(se).")
}


# ------------------------------------------------------------------------------
# Calculate uncertainty metrics
# ------------------------------------------------------------------------------
# For each TE copy in each sample, we summarize how much the
# 200 Gibbs replicates vary. High variation = high uncertainty.
#
# Four metrics:
#
# CV (Coefficient of Variation) = SD / mean
#   -> relative measure; good for comparing TEs with different
#      expression levels. A CV of 0.5 means the SD is 50% of
#      the mean.
#
# SD (Standard Deviation)
#   -> absolute spread of the replicates around the mean.
#
# IQR (Interquartile Range) = Q75 - Q25
#   -> spread of the middle 50% of replicates; robust to
#      extreme values.
#
# phi (quasi-Poisson overdispersion) = variance / mean
#   -> estimates variance inflation induced by read-to-transcript
#      mapping ambiguity. Equivalent to the moment estimator proposed
#      by Chen et al. (Nucleic Acids Research, 2024) and Smyth et al.
#      (NAR Genomics and Bioinformatics, 2024). rowVars() applies
#      Bessel's correction (divides by B-1).
#
# PSEUDO-COUNT NOTE:
# We add 1 to all replicate values before computing metrics.
# This prevents division by zero in CV (SD / mean) and phi
# (variance / mean) when a TE has zero counts in a sample.
# Adding 1 is a standard practice in genomics (pseudo-count)
# and has negligible effect on TEs with typical expression
# levels, where values are much larger than 1.

cat("Calculating uncertainty metrics...\n")

mat_cv   <- matrix(NA, nrow = n_te, ncol = n_samples,
                   dimnames = list(te_names, samples))
mat_sd   <- mat_cv
mat_iqr  <- mat_cv
mat_phi  <- mat_cv
mat_mean <- mat_cv   # mean of Gibbs replicates; used for expression profiles

for (s in seq_len(n_samples)) {
  
  cat("  Processing sample", s, "of", n_samples, ":", samples[s], "\n")
  
  rep_matrix <- sapply(rep_names, function(r) assay(se, r)[, s]) + 1
  
  mean_vals <- rowMeans(rep_matrix, na.rm = TRUE)
  var_vals  <- matrixStats::rowVars(rep_matrix, na.rm = TRUE)
  sd_vals   <- matrixStats::rowSds(rep_matrix, na.rm = TRUE)
  iqr_vals  <- matrixStats::rowIQRs(rep_matrix, na.rm = TRUE)
  
  mat_mean[, s] <- mean_vals
  mat_cv[, s]   <- sd_vals / mean_vals
  mat_sd[, s]   <- sd_vals
  mat_iqr[, s]  <- iqr_vals
  mat_phi[, s]  <- var_vals / mean_vals
}

cat("\nMetrics calculated!\n\n")


# ------------------------------------------------------------------------------
# Clean TE ID names in all metric matrices
# ------------------------------------------------------------------------------
# Some TE IDs contain the "|" character (e.g. "TE_001|extra"),
# which can break matrix and data frame operations. So we remove everything
# from "|" onward, keeping only the first part of the ID.
# If your IDs do not contain "|", this step has no effect.

clean_names <- sub("\\|.*", "", rownames(mat_cv))

rownames(mat_cv)   <- clean_names
rownames(mat_sd)   <- clean_names
rownames(mat_iqr)  <- clean_names
rownames(mat_phi)  <- clean_names
rownames(mat_mean) <- clean_names

te_names <- clean_names

cat("Rownames cleaned.\n")
cat("Example:", head(te_names, 2), "\n\n")

# Global mean CV per TE across all samples (used for Plot 01 and ranking)
global_cv_te <- rowMeans(mat_cv, na.rm = TRUE)


# ------------------------------------------------------------------------------
# Summarize metrics across samples - grouped by condition
# ------------------------------------------------------------------------------

conditions_vector <- as.character(colData(se)$condition)
condition_levels  <- c("control", "drought", "rewatering")

# Helper function: given a metric matrix and a condition label,
# return the row means across samples belonging to that condition
summarise_by_condition <- function(mat, cond) {
  cols <- which(conditions_vector == cond)
  if (length(cols) == 1) return(mat[, cols])
  rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
}

# Compute per-condition means for each metric
cv_by_cond  <- lapply(condition_levels, function(c) summarise_by_condition(mat_cv,  c))
sd_by_cond  <- lapply(condition_levels, function(c) summarise_by_condition(mat_sd,  c))
iqr_by_cond <- lapply(condition_levels, function(c) summarise_by_condition(mat_iqr, c))
phi_by_cond <- lapply(condition_levels, function(c) summarise_by_condition(mat_phi, c))

# Compute per-condition mean TPM from the abundance assay
tpm_mat     <- assay(se, "abundance")
tpm_by_cond <- lapply(condition_levels, function(c) {
  cols <- which(conditions_vector == c)
  if (length(cols) == 1) return(tpm_mat[, cols])
  rowMeans(tpm_mat[, cols, drop = FALSE], na.rm = TRUE)
})

# Name each list element by condition for clarity
names(cv_by_cond)  <- condition_levels
names(sd_by_cond)  <- condition_levels
names(iqr_by_cond) <- condition_levels
names(phi_by_cond) <- condition_levels
names(tpm_by_cond) <- condition_levels

# Strip vector names to prevent | from breaking data.frame construction
strip_names <- function(lst) lapply(lst, function(x) { names(x) <- NULL; x })
cv_by_cond  <- strip_names(cv_by_cond)
sd_by_cond  <- strip_names(sd_by_cond)
iqr_by_cond <- strip_names(iqr_by_cond)
phi_by_cond <- strip_names(phi_by_cond)
tpm_by_cond <- strip_names(tpm_by_cond)

# Build mapping from clean TE ID to base ID for classification join
id_map <- as.data.frame(
  list(
    te_copy    = te_names,
    te_base_id = sub("_copy.*", "", te_names)
  ),
  stringsAsFactors = FALSE
)
rownames(id_map) <- NULL

# Join with classification reference
id_map <- merge(
  id_map,
  te_class_ref[, c("te_id", "classification")],
  by.x  = "te_base_id",
  by.y  = "te_id",
  all.x = TRUE
)

id_map$classification[is.na(id_map$classification)] <- "Unknown"

id_map$family <- sub("/.*", "", id_map$classification)
id_map$superfamily <- ifelse(
  grepl("/", id_map$classification),
  sub(".*/", "", id_map$classification),
  id_map$family
)

id_map$family[id_map$family == "unknown"] <- "Unknown"
id_map$superfamily[id_map$superfamily == "unknown"] <- "Unknown"

# Build copy_metrics with one row per TE copy and one column set per condition
copy_metrics <- as.data.frame(
  list(
    te_copy             = te_names,
    # TPM per condition
    mean_tpm_control    = tpm_by_cond[["control"]],
    mean_tpm_drought    = tpm_by_cond[["drought"]],
    mean_tpm_rewatering = tpm_by_cond[["rewatering"]],
    # CV per condition (main uncertainty metric)
    mean_cv_control     = cv_by_cond[["control"]],
    mean_cv_drought     = cv_by_cond[["drought"]],
    mean_cv_rewatering  = cv_by_cond[["rewatering"]],
    # SD per condition
    mean_sd_control     = sd_by_cond[["control"]],
    mean_sd_drought     = sd_by_cond[["drought"]],
    mean_sd_rewatering  = sd_by_cond[["rewatering"]],
    # IQR per condition
    mean_iqr_control    = iqr_by_cond[["control"]],
    mean_iqr_drought    = iqr_by_cond[["drought"]],
    mean_iqr_rewatering = iqr_by_cond[["rewatering"]],
    # quasi-Poisson phi per condition (Chen et al. NAR 2024; Smyth et al. 2024)
    mean_phi_control    = phi_by_cond[["control"]],
    mean_phi_drought    = phi_by_cond[["drought"]],
    mean_phi_rewatering = phi_by_cond[["rewatering"]]
  ),
  stringsAsFactors = FALSE
)
rownames(copy_metrics) <- NULL

# Global mean CV across all conditions for ranking purposes
copy_metrics$mean_cv_global <- rowMeans(
  copy_metrics[, c("mean_cv_control", "mean_cv_drought", "mean_cv_rewatering")],
  na.rm = TRUE
)

# Global mean phi across all conditions
copy_metrics$mean_phi_global <- rowMeans(
  copy_metrics[, c("mean_phi_control", "mean_phi_drought", "mean_phi_rewatering")],
  na.rm = TRUE
)

# Merge classification
copy_metrics <- merge(
  copy_metrics,
  id_map[, c("te_copy", "family", "superfamily", "classification")],
  by    = "te_copy",
  all.x = TRUE
)

# Sort by global CV descending and add rank
copy_metrics <- copy_metrics[order(copy_metrics$mean_cv_global,
                                   decreasing = TRUE,
                                   na.last    = TRUE), ]
copy_metrics$rank <- seq_len(nrow(copy_metrics))

write.table(copy_metrics,
            file.path(OUT_DIR, "tables", "uncertainty_by_copy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Table saved: uncertainty_by_copy.tsv\n")
cat("Most uncertain TE:", copy_metrics$te_copy[1],
    "(global CV =", round(copy_metrics$mean_cv_global[1], 3), ")\n\n")


# ------------------------------------------------------------------------------
# Aggregate metrics by TE family
# ------------------------------------------------------------------------------

family_metrics <- copy_metrics %>%
  group_by(family, superfamily) %>%
  summarise(
    n_copies  = n(),
    mean_cv   = mean(mean_cv_global,  na.rm = TRUE),
    median_cv = median(mean_cv_global, na.rm = TRUE),
    mean_phi  = mean(mean_phi_global, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  arrange(desc(mean_cv))

write.table(family_metrics,
            file.path(OUT_DIR, "tables", "uncertainty_by_family.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Table saved: uncertainty_by_family.tsv\n\n")


# ------------------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------------------

THEME_BASE <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

condition_colors <- c(
  control    = "#56B4E9",
  drought    = "#E69F00",
  rewatering = "#009E73"
)

cat("Generating plots...\n")


# Plot 1: CV histogram — TEs vs protein-coding genes ---------------------------
# Educational comparison showing that TEs have systematically higher
# quantification uncertainty than genes, due to multi-mapping.
# X = global mean CV (averaged across all samples).
# Both curves are shown on the same axis for direct visual comparison.

global_cv_df <- data.frame(
  cv           = c(global_cv_te, global_cv_gene),
  feature_type = c(rep("Transposable Element", length(global_cv_te)),
                   rep("Protein-coding gene",  length(global_cv_gene)))
) %>% filter(!is.na(cv))

med_df <- global_cv_df %>%
  group_by(feature_type) %>%
  summarise(med = median(cv, na.rm = TRUE), .groups = "drop")

p01 <- ggplot(global_cv_df,
              aes(x = cv, fill = feature_type, colour = feature_type)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, alpha = 0.45, position = "identity") +
  geom_density(alpha = 0, linewidth = 0.8) +
  geom_vline(data = med_df,
             aes(xintercept = med, colour = feature_type),
             linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values   = c("Transposable Element" = "#D55E00",
                                 "Protein-coding gene"  = "#0072B2")) +
  scale_colour_manual(values = c("Transposable Element" = "#D55E00",
                                 "Protein-coding gene"  = "#0072B2")) +
  labs(
    title    = "Plot 01: Quantification uncertainty — TEs vs Protein-coding Genes",
    subtitle = paste0("Mean CV across all ", n_samples, " samples; ",
                      n_te, " TEs and ", n_gene, " genes; ",
                      "dashed lines = medians"),
    x        = "Mean CV  (SD / mean of Gibbs replicates)",
    y        = "Density",
    fill     = "Feature type",
    colour   = "Feature type"
  ) +
  THEME_BASE

ggsave(file.path(OUT_DIR, "plots", "01_cv_histogram_TE_vs_gene.pdf"),
       p01, width = 10, height = 5)
cat("  -> 01_cv_histogram_TE_vs_gene.pdf saved\n")


# Plot 2: CV histogram by condition (TEs only) ---------------------------------
# Shows whether quantification uncertainty shifts between conditions.
# Each curve aggregates all TE copies from samples in that condition.
# Dashed lines mark the median per condition.

cv_cond_long <- data.frame(
  condition = rep(condition_levels, each = nrow(copy_metrics)),
  cv        = c(copy_metrics$mean_cv_control,
                copy_metrics$mean_cv_drought,
                copy_metrics$mean_cv_rewatering)
) %>%
  filter(!is.na(cv)) %>%
  mutate(condition = factor(condition, levels = condition_levels))

med_cond <- cv_cond_long %>%
  group_by(condition) %>%
  summarise(med = median(cv, na.rm = TRUE), .groups = "drop")

p02 <- ggplot(cv_cond_long,
              aes(x = cv, fill = condition, colour = condition)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 80, alpha = 0.4, position = "identity") +
  geom_density(alpha = 0, linewidth = 0.8) +
  geom_vline(data = med_cond,
             aes(xintercept = med, colour = condition),
             linetype = "dashed", linewidth = 0.9) +
  scale_fill_manual(values   = condition_colors) +
  scale_colour_manual(values = condition_colors) +
  labs(
    title    = "Plot 02: TE quantification uncertainty by condition",
    subtitle = paste0(n_te, " TE copies; dashed lines = medians per condition"),
    x        = "Mean CV",
    y        = "Density",
    fill     = "Condition",
    colour   = "Condition"
  ) +
  THEME_BASE

ggsave(file.path(OUT_DIR, "plots", "02_cv_histogram_by_condition.pdf"),
       p02, width = 9, height = 5)
cat("  -> 02_cv_histogram_by_condition.pdf saved\n")


# Plot 3: CV violin + boxplot by condition (TEs only) --------------------------
# Directly compares uncertainty distributions across the three conditions.
# The violin shows the full distribution; the embedded boxplot shows
# median and IQR.

p03 <- ggplot(cv_cond_long,
              aes(x = condition, y = cv, fill = condition)) +
  geom_violin(trim = TRUE, alpha = 0.5, scale = "width") +
  geom_boxplot(width = 0.15, fill = "white",
               outlier.size = 0.4, outlier.alpha = 0.3) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title    = "Plot 03: CV distribution per condition (TEs)",
    subtitle = "Violin + boxplot aggregating all TE copies",
    x        = "Condition",
    y        = "Mean CV",
    fill     = "Condition"
  ) +
  THEME_BASE +
  theme(legend.position = "none")

ggsave(file.path(OUT_DIR, "plots", "03_cv_boxplot_by_condition.pdf"),
       p03, width = 7, height = 5)
cat("  -> 03_cv_boxplot_by_condition.pdf saved\n")


# Plot 4: Expression profiles with error bars — top N uncertain TEs ------------
# For each of the top N most uncertain TEs, one panel shows the expression
# profile across all samples, with error bars representing the uncertainty.
#
# X  = sample (grouped and coloured by condition)
# Y  = mean of scaled Gibbs replicates (expression estimate)
# Error bars = ±1 SD of Gibbs replicates (quantification uncertainty)
#
# A wide error bar means the Salmon posterior is spread — the true count
# for that TE in that sample is poorly constrained.

top_te_ids <- head(copy_metrics$te_copy, N_TOP_TE_PROFILES)

profile_list <- list()

for (te in top_te_ids) {
  idx <- which(te_names == te)
  if (length(idx) == 0) next
  profile_list[[te]] <- data.frame(
    te_copy   = te,
    sample    = samples,
    condition = conditions_vector,
    mean_expr = as.numeric(mat_mean[idx, ]),
    sd_expr   = as.numeric(mat_sd[idx, ]),
    stringsAsFactors = FALSE
  )
}

profile_df <- bind_rows(profile_list) %>%
  mutate(
    condition = factor(condition, levels = condition_levels),
    sample    = factor(sample,
                       levels = samples[order(conditions_vector,
                                              method = "radix")])
  )

p05 <- ggplot(profile_df,
              aes(x = sample, y = mean_expr,
                  colour = condition, group = te_copy)) +
  geom_errorbar(aes(ymin = mean_expr - sd_expr,
                    ymax = mean_expr + sd_expr),
                width = 0.3, alpha = 0.7) +
  geom_point(size = 2.5) +
  geom_line(alpha = 0.5) +
  scale_colour_manual(values = condition_colors) +
  facet_wrap(~ te_copy, scales = "free_y",
             ncol = 3, labeller = label_wrap_gen(25)) +
  labs(
    title    = paste0("Plot 04: Expression profiles — top ",
                      N_TOP_TE_PROFILES, " most uncertain TEs"),
    subtitle = "Points = mean of Gibbs replicates; error bars = \u00b11 SD",
    x        = "Sample",
    y        = "Mean expression (scaled Gibbs replicates)",
    colour   = "Condition"
  ) +
  THEME_BASE +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text  = element_text(size = 7))

ggsave(file.path(OUT_DIR, "plots", "04_expression_profiles_top_TEs.pdf"),
       p05,
       width  = 14,
       height = ceiling(N_TOP_TE_PROFILES / 3) * 3.5 + 2)
cat("  -> 04_expression_profiles_top_TEs.pdf saved\n\n")


# ------------------------------------------------------------------------------
# Final summary
# ------------------------------------------------------------------------------

n_complete_te <- sum(complete.cases(mat_cv))

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("TEs analysed:              ", n_te,      "\n")
cat("Protein-coding genes:      ", n_gene,    "\n")
cat("Samples:                   ", n_samples, "\n")
cat("Gibbs replicates:          ", n_gibbs,   "\n")
cat("Families found:            ", length(unique(copy_metrics$family)), "\n")
cat("TEs with complete CV data: ", n_complete_te,
    "(", round(n_complete_te / n_te * 100, 1), "% of total)\n")
cat("\nMean CV by condition:\n")
for (cond in condition_levels) {
  col <- paste0("mean_cv_", cond)
  cat("  ", formatC(cond, width = 12, flag = "-"), ":",
      round(mean(copy_metrics[[col]], na.rm = TRUE), 4), "\n")
}
cat("\nGlobal mean CV — TEs:  ",
    round(mean(global_cv_te,   na.rm = TRUE), 4), "\n")
cat("Global mean CV — genes:",
    round(mean(global_cv_gene, na.rm = TRUE), 4), "\n")
cat("\nMost uncertain TE:         ", copy_metrics$te_copy[1],
    "(CV =", round(copy_metrics$mean_cv_global[1], 4), ")\n")
cat("\nOutput files saved in:", OUT_DIR, "\n")
cat("  tables/\n")
cat("    uncertainty_by_copy.tsv\n")
cat("    uncertainty_by_family.tsv\n")
cat("  plots/\n")
cat("    01_cv_histogram_TE_vs_gene.pdf\n")
cat("    02_cv_histogram_by_condition.pdf\n")
cat("    03_cv_boxplot_by_condition.pdf\n")
cat("    04_expression_profiles_top_TEs.pdf\n")
cat("=======================================================================\n")

sessionInfo()
