# ==============================================================================
# TE quantification uncertainty analysis using Salmon Gibbs replicates.
# VERSION: single-sample mode (for testing)
#
# The infReps structure
# is correctly handled for both single and multi-sample runs.
#
# Usage (terminal):
# Rscript quantification_uncertainty_single_sample.R <conditions.tsv> [salmon_dir] [out_dir]
# ==============================================================================


# ==============================================================================
# STEP 1: Install and load packages
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("ggplot2", "dplyr", "tidyr", "matrixStats",
                   "pheatmap", "viridis")

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
library(pheatmap)
library(viridis)


# ==============================================================================
# STEP 2: Set paths
# ==============================================================================

CONDITIONS_FILE <- "~/rnaseq_d.tsv"
SALMON_DIR      <- "~/salmon_results"
OUT_DIR         <- "~/test_results"

if (!file.exists(CONDITIONS_FILE)) {
  stop("Conditions file not found: ", CONDITIONS_FILE)
}


# ==============================================================================
# STEP 3: Read the conditions file
# ==============================================================================

conditions <- read.table(
  CONDITIONS_FILE,
  header           = TRUE,
  sep              = "\t",
  stringsAsFactors = FALSE
)

# --- Diagnostic: show what was loaded ---
cat("Conditions file loaded successfully!\n")
cat("Number of samples:", nrow(conditions), "\n")
cat("Columns found:", paste(colnames(conditions), collapse = ", "), "\n")
cat("Groups found:", paste(unique(conditions$condition), collapse = ", "), "\n\n")

# Single-sample mode warning
if (nrow(conditions) == 1) {
  cat("*** SINGLE-SAMPLE MODE: some multi-sample plots will be skipped. ***\n\n")
}

dir.create(file.path(OUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "plots"),  recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# STEP 4: Locate Salmon output files
# ==============================================================================

conditions$files <- file.path(SALMON_DIR, conditions$sample, "quant.sf")

missing_files <- !file.exists(conditions$files)
if (any(missing_files)) {
  stop("The following quant.sf files were not found:\n",
       paste(conditions$files[missing_files], collapse = "\n"))
}

cat("All quant.sf files found!\n")
cat("Files to be imported:\n")
print(conditions$files)
cat("\n")

files <- setNames(conditions$files, conditions$sample)


# ==============================================================================
# STEP 5: Import Salmon data with tximport
# ==============================================================================

cat("Importing Salmon data (this may take a few minutes)...\n")

txi <- tximport(
  files,
  type        = "salmon",
  txOut       = TRUE,
  dropInfReps = FALSE
)

# --- Diagnostic: check dimensions before building SE ---
cat("\n--- Dimension check ---\n")
cat("counts:   ", dim(txi$counts), "\n")
cat("abundance:", dim(txi$abundance), "\n")
cat("length:   ", dim(txi$length), "\n")
cat("infReps is a list of", length(txi$infReps), "element(s)\n")
cat("Each infRep element dimensions:", dim(txi$infReps[[1]]), "\n")
cat("  (expected: n_TEs x n_GibbsSamples)\n\n")

n_te_raw  <- nrow(txi$counts)
n_samples <- ncol(txi$counts)
n_gibbs   <- ncol(txi$infReps[[1]])   # number of Gibbs replicates (e.g. 200)

cat("TEs imported:          ", n_te_raw, "\n")
cat("Samples detected:      ", n_samples, "\n")
cat("Gibbs replicates/TE:   ", n_gibbs, "\n\n")

if (n_gibbs == 0) {
  stop("No Gibbs replicates found. Make sure Salmon was run with --numGibbsSamples.")
}


# ==============================================================================
# STEP 6: Build a SummarizedExperiment object
# ==============================================================================
# tximport returns txi$infReps as a LIST with one element per SAMPLE.
# Each element is a matrix of n_TEs x n_GibbsSamples.
#
# SummarizedExperiment needs each assay to be n_TEs x n_SAMPLES.
# So we restructure: for Gibbs replicate k, collect column k from
# every sample and cbind them into one n_TEs x n_SAMPLES matrix.
# ==============================================================================

cat("Restructuring infReps for SummarizedExperiment...\n")

inf_assays <- setNames(
  lapply(seq_len(n_gibbs), function(k) {
    # For replicate k: extract column k from each sample's matrix and bind
    mat <- do.call(cbind, lapply(txi$infReps, function(sample_mat) {
      sample_mat[, k, drop = FALSE]
    }))
    colnames(mat) <- conditions$sample
    mat
  }),
  paste0("infRep", seq_len(n_gibbs))
)

# Verify the restructuring worked
cat("infRep1 dimensions after restructuring:", dim(inf_assays[[1]]), "\n")
cat("(expected:", n_te_raw, "x", n_samples, ")\n\n")

# Build the SummarizedExperiment
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

# Scale inferential replicates (adjusts for sequencing depth)
# scaleInfReps normalizes across samples; skip when only one sample is present
if (n_samples > 1) {
  se <- scaleInfReps(se)
} else {
  cat("scaleInfReps skipped: normalization across samples requires at least 2 samples.\n\n")
}

samples   <- colnames(se)
n_samples <- length(samples)
te_names  <- rownames(se)
n_te      <- nrow(se)

cat("SummarizedExperiment built successfully!\n")
cat("Dimensions:", n_te, "TEs x", n_samples, "samples\n\n")


# ==============================================================================
# STEP 7: Extract TE family and superfamily names
# ==============================================================================
# TE names from RepeatMasker: CopyName#Family/Superfamily
# Example: "AluSx1#SINE/Alu"
# ==============================================================================

te_family <- sub(".*#", "", te_names)
te_family <- sub("/.*", "", te_family)
te_family[!grepl("#", te_names)] <- "unknown"

te_superfamily <- rep("unknown", length(te_names))
has_hash  <- grepl("#", te_names)
has_slash <- grepl("/", te_names)
te_superfamily[has_hash & has_slash] <- sub(".*#.*?/", "", te_names[has_hash & has_slash])

cat("Unique families found:     ", length(unique(te_family)), "\n")
cat("Unique superfamilies found:", length(unique(te_superfamily)), "\n\n")


# ==============================================================================
# STEP 8: Calculate uncertainty metrics
# ==============================================================================
# For each TE in each sample we summarize variance across Gibbs replicates.
#
# CV   = SD / mean        (relative spread; main metric)
# SD   = standard deviation
# IQR  = Q75 - Q25        (robust to outliers)
# Fano = variance / mean  (suited for count data)
# ==============================================================================

cat("Calculating uncertainty metrics...\n")

rep_names <- paste0("infRep", seq_len(n_gibbs))

mat_cv   <- matrix(NA, nrow = n_te, ncol = n_samples,
                   dimnames = list(te_names, samples))
mat_sd   <- mat_cv
mat_iqr  <- mat_cv
mat_fano <- mat_cv

for (s in seq_len(n_samples)) {
  
  cat("  Processing sample", s, "of", n_samples, ":", samples[s], "\n")
  
  # Extract all Gibbs replicates for this sample -> matrix n_TEs x n_gibbs
  rep_matrix <- sapply(rep_names, function(r) assay(se, r)[, s])
  
  mean_vals <- rowMeans(rep_matrix, na.rm = TRUE)
  var_vals  <- matrixStats::rowVars(rep_matrix, na.rm = TRUE)
  sd_vals   <- matrixStats::rowSds(rep_matrix, na.rm = TRUE)
  iqr_vals  <- matrixStats::rowIQRs(rep_matrix, na.rm = TRUE)
  
  mat_cv[, s]   <- ifelse(mean_vals > 0, sd_vals / mean_vals, NA)
  mat_sd[, s]   <- sd_vals
  mat_iqr[, s]  <- iqr_vals
  mat_fano[, s] <- ifelse(mean_vals > 0, var_vals / mean_vals, NA)
}

cat("\nMetrics calculated!\n\n")


# ==============================================================================
# STEP 9: Summarize metrics across samples
# ==============================================================================

mean_cv   <- rowMeans(mat_cv,   na.rm = TRUE)
mean_sd   <- rowMeans(mat_sd,   na.rm = TRUE)
mean_iqr  <- rowMeans(mat_iqr,  na.rm = TRUE)
mean_fano <- rowMeans(mat_fano, na.rm = TRUE)
mean_tpm  <- rowMeans(assay(se, "abundance"), na.rm = TRUE)

copy_metrics <- data.frame(
  te_copy     = te_names,
  family      = te_family,
  superfamily = te_superfamily,
  mean_tpm    = mean_tpm,
  mean_cv     = mean_cv,
  mean_sd     = mean_sd,
  mean_iqr    = mean_iqr,
  mean_fano   = mean_fano,
  stringsAsFactors = FALSE
)

copy_metrics <- copy_metrics[order(copy_metrics$mean_cv,
                                   decreasing = TRUE,
                                   na.last    = TRUE), ]
copy_metrics$rank <- seq_len(nrow(copy_metrics))

write.table(copy_metrics,
            file.path(OUT_DIR, "tables", "uncertainty_by_copy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Table saved: uncertainty_by_copy.tsv\n")
cat("Most uncertain TE:", copy_metrics$te_copy[1],
    "(CV =", round(copy_metrics$mean_cv[1], 3), ")\n\n")


# ==============================================================================
# STEP 10: Aggregate metrics by TE family
# ==============================================================================

family_metrics <- copy_metrics %>%
  group_by(family, superfamily) %>%
  summarise(
    n_copies   = n(),
    median_tpm = median(mean_tpm, na.rm = TRUE),
    mean_cv    = mean(mean_cv,    na.rm = TRUE),
    median_cv  = median(mean_cv,  na.rm = TRUE),
    mean_fano  = mean(mean_fano,  na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  arrange(desc(mean_cv))

write.table(family_metrics,
            file.path(OUT_DIR, "tables", "uncertainty_by_family.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Table saved: uncertainty_by_family.tsv\n\n")


# ==============================================================================
# STEP 11: Plots
# ==============================================================================

THEME_BASE <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

cat("Generating plots...\n")


# Plot 1: CV distribution (density) -------------------------------------------

p1 <- ggplot(copy_metrics, aes(x = mean_cv)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(
    xintercept = median(copy_metrics$mean_cv, na.rm = TRUE),
    linetype   = "dashed", colour = "firebrick", linewidth = 0.9
  ) +
  annotate(
    "text",
    x      = median(copy_metrics$mean_cv, na.rm = TRUE) * 1.05,
    y      = Inf, vjust = 2, hjust = 0, size = 3.5,
    label  = paste0("Median = ", round(median(copy_metrics$mean_cv, na.rm = TRUE), 3)),
    colour = "firebrick"
  ) +
  labs(
    title    = "Plot 1: Global CV distribution",
    subtitle = paste0(n_te, " TE copies — ", n_samples, " sample(s)"),
    x        = "Mean CV  (SD / mean of Gibbs replicates)",
    y        = "Density"
  ) +
  THEME_BASE

ggsave(file.path(OUT_DIR, "plots", "01_cv_distribution.pdf"),
       p1, width = 8, height = 5)
cat("  -> 01_cv_distribution.pdf saved\n")


# Plot 2: Expression level vs CV ----------------------------------------------

p2 <- ggplot(copy_metrics,
             aes(x = log10(mean_tpm + 0.01), y = mean_cv, colour = superfamily)) +
  geom_point(size = 0.8, alpha = 0.4) +
  geom_smooth(method = "loess", colour = "black", se = TRUE, linewidth = 0.8) +
  labs(
    title    = "Plot 2: Expression level vs Uncertainty",
    subtitle = "Each dot = one TE copy; black line = LOESS trend",
    x        = "log10(mean TPM + 0.01)",
    y        = "Mean CV",
    colour   = "Superfamily"
  ) +
  THEME_BASE +
  theme(legend.text = element_text(size = 7))

ggsave(file.path(OUT_DIR, "plots", "02_abundance_vs_cv.pdf"),
       p2, width = 9, height = 6)
cat("  -> 02_abundance_vs_cv.pdf saved\n")


# Plot 3: Boxplot — CV by family (top 25) -------------------------------------

top25_families <- family_metrics %>%
  arrange(desc(n_copies)) %>%
  slice_head(n = 25) %>%
  pull(family)

copy_top25 <- copy_metrics %>% filter(family %in% top25_families)

p3 <- ggplot(copy_top25,
             aes(x = reorder(family, mean_cv, FUN = median), y = mean_cv)) +
  geom_boxplot(fill = "steelblue", alpha = 0.6,
               outlier.size = 0.5, outlier.alpha = 0.3) +
  coord_flip() +
  labs(
    title    = "Plot 3: CV distribution by TE family",
    subtitle = "Top 25 families by number of copies; ordered by median CV",
    x        = "Family",
    y        = "Mean CV"
  ) +
  THEME_BASE

ggsave(file.path(OUT_DIR, "plots", "03_cv_by_family_boxplot.pdf"),
       p3, width = 9, height = 9)
cat("  -> 03_cv_by_family_boxplot.pdf saved\n")


# Plot 4: Heatmap — CV per sample, top 40 TEs ---------------------------------

top40_te     <- head(copy_metrics$te_copy, 40)
mat_cv_top40 <- mat_cv[top40_te, , drop = FALSE]

# Remove rows with any NA to allow hclust to run
mat_cv_top40 <- mat_cv_top40[complete.cases(mat_cv_top40), , drop = FALSE]

# If too few rows remain, disable row clustering to avoid hclust error
cluster_rows_flag <- nrow(mat_cv_top40) > 1

annotation_columns <- data.frame(
  Condition = as.character(colData(se)$condition),
  row.names = colnames(se)
)

pdf(file.path(OUT_DIR, "plots", "04_heatmap_cv_top40_TEs.pdf"),
    width = max(8, n_samples * 0.7 + 3), height = 10)

pheatmap(
  mat_cv_top40,
  color          = viridis(100),
  cluster_rows   = cluster_rows_flag,
  cluster_cols   = (n_samples > 1),
  annotation_col = annotation_columns,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  fontsize_row   = 7,
  main           = "CV per sample — Top 40 most uncertain TEs",
  na_col         = "grey80"
)

dev.off()
cat("  -> 04_heatmap_cv_top40_TEs.pdf saved\n")

# Plot 5: Violin — CV per sample (skipped if only one sample) -----------------

if (n_samples > 1) {
  
  cv_long <- as.data.frame(mat_cv) %>%
    mutate(te_copy = te_names) %>%
    pivot_longer(
      cols      = all_of(samples),
      names_to  = "sample",
      values_to = "cv"
    ) %>%
    left_join(conditions[, c("sample", "condition")], by = "sample") %>%
    filter(!is.na(cv))
  
  p5 <- ggplot(cv_long, aes(x = sample, y = cv, fill = condition)) +
    geom_violin(trim = TRUE, alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white",
                 outlier.size = 0.3, alpha = 0.9) +
    labs(
      title    = "Plot 5: CV distribution per sample",
      subtitle = "Each violin = all TEs in that sample",
      x        = "Sample",
      y        = "CV (SD / mean of Gibbs replicates)",
      fill     = "Condition"
    ) +
    THEME_BASE +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUT_DIR, "plots", "05_cv_violin_per_sample.pdf"),
         p5, width = max(8, n_samples * 0.8 + 2), height = 6)
  cat("  -> 05_cv_violin_per_sample.pdf saved\n\n")
  
} else {
  cat("  -> Plot 5 (violin per sample) skipped: only one sample available.\n\n")
}


# ==============================================================================
# STEP 12: Final summary
# ==============================================================================

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("TEs analysed:        ", n_te, "\n")
cat("Samples:             ", n_samples, "\n")
cat("Gibbs replicates:    ", n_gibbs, "\n")
cat("Families found:      ", length(unique(te_family)), "\n")
cat("\nGlobal mean CV:      ", round(mean(copy_metrics$mean_cv, na.rm = TRUE), 4), "\n")
cat("Global median CV:    ", round(median(copy_metrics$mean_cv, na.rm = TRUE), 4), "\n")
cat("Most uncertain TE:   ", copy_metrics$te_copy[1],
    "(CV =", round(copy_metrics$mean_cv[1], 4), ")\n")
cat("\nOutput files saved in:", OUT_DIR, "\n")
cat("  tables/\n")
cat("    uncertainty_by_copy.tsv\n")
cat("    uncertainty_by_family.tsv\n")
cat("  plots/\n")
cat("    01_cv_distribution.pdf\n")
cat("    02_abundance_vs_cv.pdf\n")
cat("    03_cv_by_family_boxplot.pdf\n")
cat("    04_heatmap_cv_top40_TEs.pdf\n")
if (n_samples > 1) cat("    05_cv_violin_per_sample.pdf\n")
cat("=======================================================================\n")