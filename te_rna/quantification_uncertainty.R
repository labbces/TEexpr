# ==============================================================================
# TE quantification uncertainty analysis using Salmon Gibbs replicates.
#
# Compatible with single and multiple samples.
#
# Usage (terminal):
# Rscript quantification_uncertainty_final.R
# ==============================================================================


# ==============================================================================
# STEP 1: Install and load packages
# ==============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c("ggplot2", "dplyr", "tidyr", "matrixStats",
                   "pheatmap", "viridis", "gridExtra", "gtable", "grid")

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
library(gridExtra)
library(gtable)
library(grid)


# ==============================================================================
# STEP 2: Set paths
# ==============================================================================

CONDITIONS_FILE     <- "~/rnaseq_d.tsv"
SALMON_DIR          <- "~/salmon_results"
OUT_DIR             <- "~/uncertainty_results"
CLASSIFICATION_FILE <- "~/Sviridis.flTE.mapids"

if (!file.exists(CONDITIONS_FILE)) {
  stop("Conditions file not found: ", CONDITIONS_FILE)
}

if (!file.exists(CLASSIFICATION_FILE)) {
  stop("Classification file not found: ", CLASSIFICATION_FILE)
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

cat("Conditions file loaded!\n")
cat("Number of samples:", nrow(conditions), "\n")
cat("Columns found:", paste(colnames(conditions), collapse = ", "), "\n")
cat("Groups found:", paste(unique(conditions$condition), collapse = ", "), "\n\n")

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

cat("All quant.sf files found!\n\n")

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

cat("\n--- Dimension check after import ---\n")
cat("counts:   ", dim(txi$counts),    "\n")
cat("abundance:", dim(txi$abundance), "\n")
cat("length:   ", dim(txi$length),    "\n")
cat("infReps: list of", length(txi$infReps), "element(s)\n")
cat("Each infRep element:", dim(txi$infReps[[1]]), "(TEs x Gibbs replicates)\n\n")

n_gibbs <- ncol(txi$infReps[[1]])

cat("TEs imported:         ", nrow(txi$counts), "\n")
cat("Samples detected:     ", ncol(txi$counts), "\n")
cat("Gibbs replicates/TE:  ", n_gibbs, "\n\n")

if (n_gibbs == 0) {
  stop("No Gibbs replicates found. Make sure Salmon was run with --numGibbsSamples.")
}


# ==============================================================================
# STEP 6: Build a SummarizedExperiment object
# ==============================================================================
# tximport returns txi$infReps as a list indexed by sample.
# Each element is a matrix n_TEs x n_Gibbs.
# SE needs each assay as n_TEs x n_samples.
# We restructure: for replicate k, collect column k from every sample.
# ==============================================================================

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

# Scale inferential replicates (requires at least 2 samples)
if (ncol(se) > 1) {
  se <- scaleInfReps(se)
} else {
  cat("scaleInfReps skipped: requires at least 2 samples.\n\n")
}

samples   <- colnames(se)
n_samples <- length(samples)

cat("SummarizedExperiment built!\n")
cat("Dimensions:", nrow(se), "IDs x", n_samples, "samples\n\n")


# ==============================================================================
# STEP 7: Assign TE classification from external file
# ==============================================================================
# TE IDs follow: TE_<numericID>_copy<N>|Chr_XX:start-end|strand
# Genes are IDs NOT starting with "TE_".
# TEdistill "[specie].mapids" file: no header, tab-separated, 3 columns:
#   col 1: TE_<numericID>
#   col 2: ClassificationA/a  <- used here
#   col 3: TE_<numericID>#ClassificationA/a
# ==============================================================================

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

# Verify extraction
cat("Extraction check:\n")
cat("  Input  :", head(te_ids, 2),      "\n")
cat("  Extracted:", head(te_base_id, 2), "\n\n")

# Check coverage
n_matched   <- sum(te_base_id %in% te_class_ref$te_id)
n_unmatched <- sum(!te_base_id %in% te_class_ref$te_id)
cat("Matched:  ", n_matched,   "\n")
cat("Unmatched:", n_unmatched, "\n\n")

# Map classification using lookup
class_lookup     <- setNames(te_class_ref$classification, te_class_ref$te_id)
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

cat("Unique classifications:", length(unique(te_classification)), "\n")
cat("Unique families:       ", length(unique(te_family)),         "\n")
cat("Unique superfamilies:  ", length(unique(te_superfamily)),    "\n\n")


# ==============================================================================
# STEP 7b: Filter SE to retain only TE rows
# ==============================================================================

se       <- se[is_te, ]
te_names <- rownames(se)
n_te     <- nrow(se)

cat("SE filtered to TEs only.\n")
cat("Dimensions:", n_te, "TEs x", n_samples, "samples\n\n")

# Alignment check
cat("Alignment check:\n")
cat("  length(te_family):     ", length(te_family),     "\n")
cat("  length(te_superfamily):", length(te_superfamily), "\n")
cat("  nrow(se):              ", nrow(se),              "\n\n")

if (length(te_family) != nrow(se) || length(te_superfamily) != nrow(se)) {
  stop("Alignment error: te_family or te_superfamily length does not match nrow(se).")
}


# ==============================================================================
# STEP 8: Calculate uncertainty metrics
# ==============================================================================
# For each TE in each sample, summarize variance across 200 Gibbs replicates.
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
# STEP 8b: Sanitize rownames of all metric matrices
# ------------------------------------------------------------------------------
# The | character in TE IDs breaks R matrix and data.frame operations.
# Truncate all rownames here once so all downstream steps work cleanly.
# te_names is also updated to stay consistent.
# ==============================================================================

clean_names <- sub("\\|.*", "", rownames(mat_cv))

rownames(mat_cv)   <- clean_names
rownames(mat_sd)   <- clean_names
rownames(mat_iqr)  <- clean_names
rownames(mat_fano) <- clean_names

te_names <- clean_names

cat("Rownames sanitized.\n")
cat("Example:", head(te_names, 2), "\n\n")


# ==============================================================================
# STEP 9: Summarize metrics across samples
# ==============================================================================

mean_cv   <- rowMeans(mat_cv,   na.rm = TRUE)
mean_sd   <- rowMeans(mat_sd,   na.rm = TRUE)
mean_iqr  <- rowMeans(mat_iqr,  na.rm = TRUE)
mean_fano <- rowMeans(mat_fano, na.rm = TRUE)
mean_tpm  <- rowMeans(assay(se, "abundance"), na.rm = TRUE)

# Strip names to prevent | from breaking data.frame construction
names(mean_cv)   <- NULL
names(mean_sd)   <- NULL
names(mean_iqr)  <- NULL
names(mean_fano) <- NULL
names(mean_tpm)  <- NULL

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

id_map$family      <- sub("/.*", "", id_map$classification)
id_map$superfamily <- ifelse(
  grepl("/", id_map$classification),
  sub(".*/", "", id_map$classification),
  id_map$family
)

id_map$family[id_map$family           == "unknown"] <- "Unknown"
id_map$superfamily[id_map$superfamily == "unknown"] <- "Unknown"

# Build copy_metrics via list to avoid | parsing issues
copy_metrics <- as.data.frame(
  list(
    te_copy   = te_names,
    mean_tpm  = mean_tpm,
    mean_cv   = mean_cv,
    mean_sd   = mean_sd,
    mean_iqr  = mean_iqr,
    mean_fano = mean_fano
  ),
  stringsAsFactors = FALSE
)
rownames(copy_metrics) <- NULL

# Merge classification
copy_metrics <- merge(
  copy_metrics,
  id_map[, c("te_copy", "family", "superfamily", "classification")],
  by    = "te_copy",
  all.x = TRUE
)

# Sort by CV descending and add rank
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
# STEP 9b: Hybrid uncertainty metric
# ------------------------------------------------------------------------------
# CV is unstable when mean TPM is near zero.
# Below TPM_THRESHOLD we use Fano index instead.
# ==============================================================================

TPM_THRESHOLD <- 0.5

copy_metrics$metric_used <- ifelse(
  copy_metrics$mean_tpm >= TPM_THRESHOLD, "CV", "Fano"
)

copy_metrics$hybrid_uncertainty <- ifelse(
  copy_metrics$metric_used == "CV",
  copy_metrics$mean_cv,
  copy_metrics$mean_fano
)

n_cv_regime   <- sum(copy_metrics$metric_used == "CV",   na.rm = TRUE)
n_fano_regime <- sum(copy_metrics$metric_used == "Fano", na.rm = TRUE)

cat("TEs using CV   (mean TPM >=", TPM_THRESHOLD, "):", n_cv_regime,   "\n")
cat("TEs using Fano (mean TPM < ", TPM_THRESHOLD, "):", n_fano_regime, "\n\n")


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


# Plot 1: Hybrid uncertainty distribution (two panels) ------------------------

p1a <- ggplot(
  filter(copy_metrics, metric_used == "CV"),
  aes(x = mean_cv)
) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  geom_vline(
    xintercept = median(copy_metrics$mean_cv[copy_metrics$metric_used == "CV"],
                        na.rm = TRUE),
    linetype = "dashed", colour = "firebrick", linewidth = 0.9
  ) +
  labs(
    title    = "1A: CV distribution (expressed TEs)",
    subtitle = paste0("mean TPM >= ", TPM_THRESHOLD, " — n = ", n_cv_regime, " TEs"),
    x = "Mean CV",
    y = "Density"
  ) +
  THEME_BASE

p1b <- ggplot(
  filter(copy_metrics, metric_used == "Fano"),
  aes(x = mean_fano)
) +
  geom_density(fill = "darkorange", alpha = 0.5) +
  geom_vline(
    xintercept = median(copy_metrics$mean_fano[copy_metrics$metric_used == "Fano"],
                        na.rm = TRUE),
    linetype = "dashed", colour = "firebrick", linewidth = 0.9
  ) +
  labs(
    title    = "1B: Fano distribution (lowly expressed TEs)",
    subtitle = paste0("mean TPM < ", TPM_THRESHOLD, " — n = ", n_fano_regime, " TEs"),
    x = "Mean Fano index",
    y = "Density"
  ) +
  THEME_BASE

pdf(file.path(OUT_DIR, "plots", "01_uncertainty_distribution.pdf"),
    width = 8, height = 9)
gridExtra::grid.arrange(p1a, p1b, ncol = 1)
dev.off()
cat("  -> 01_uncertainty_distribution.pdf saved\n")


# Plot 2: Expression level vs hybrid uncertainty ------------------------------

top_sf <- names(sort(table(copy_metrics$superfamily), decreasing = TRUE))[1:8]

copy_metrics$superfamily_plot <- ifelse(
  copy_metrics$superfamily %in% top_sf,
  copy_metrics$superfamily,
  "Other"
)

sf_order <- c(sort(top_sf), "Other")
copy_metrics$superfamily_plot <- factor(copy_metrics$superfamily_plot,
                                        levels = sf_order)

sf_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#444444"
)
names(sf_colors) <- sf_order

p2 <- ggplot(
  copy_metrics,
  aes(x = log10(mean_tpm + 0.01), y = hybrid_uncertainty,
      colour = superfamily_plot)
) +
  geom_point(size = 0.8, alpha = 0.35) +
  geom_smooth(method = "loess", colour = "black", se = TRUE, linewidth = 0.8) +
  geom_vline(
    xintercept = log10(TPM_THRESHOLD + 0.01),
    linetype = "dashed", colour = "grey30", linewidth = 0.8
  ) +
  annotate(
    "text",
    x = log10(TPM_THRESHOLD + 0.01) + 0.05,
    y = Inf, vjust = 2, hjust = 0, size = 3.2, colour = "grey30",
    label = paste0("Threshold = ", TPM_THRESHOLD, " TPM")
  ) +
  scale_colour_manual(values = sf_colors) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(
    title    = "Plot 2: Expression level vs Uncertainty",
    subtitle = "Top 8 superfamilies shown; remaining collapsed into 'Other'",
    x        = "log10(mean TPM + 0.01)",
    y        = "Hybrid uncertainty metric",
    colour   = "Superfamily"
  ) +
  THEME_BASE +
  theme(legend.text = element_text(size = 9))

ggsave(file.path(OUT_DIR, "plots", "02_abundance_vs_uncertainty.pdf"),
       p2, width = 10, height = 6)
cat("  -> 02_abundance_vs_uncertainty.pdf saved\n")


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


# Plot 4: Heatmap — CV per sample (top 40 per sample) -------------------------

# ==============================================================================
# CUSTOMIZE APPEARANCE HERE
# ==============================================================================

condition_colors <- list(
  " " = c(
    "control" = "#56B4E9",
    "drought" = "#E69F00"
  )
)

title_fontsize     <- 13        # font size in points
title_fontface     <- "bold"    # "plain", "bold", "italic", "bold.italic"
title_margin_lines <- 2         # blank lines between title and heatmap body
cell_border_color  <- "white"   # "NA" = no border

# ==============================================================================

for (s in samples) {
  
  cv_sample <- mat_cv[, s]
  valid_te  <- names(cv_sample)[!is.na(cv_sample) & !is.nan(cv_sample) &
                                  !is.infinite(cv_sample)]
  
  cat("Sample", s, "— TEs with valid CV:", length(valid_te), "\n")
  
  if (length(valid_te) < 2) {
    cat("  -> Heatmap skipped for", s, ": fewer than 2 valid TEs.\n")
    next
  }
  
  top_n    <- min(40, length(valid_te))
  top40_te <- names(sort(cv_sample[valid_te], decreasing = TRUE))[1:top_n]
  
  mat_sub <- mat_cv[top40_te, , drop = FALSE]
  
  # Impute remaining NAs with row mean
  mat_sub <- t(apply(mat_sub, 1, function(row) {
    row[is.na(row)] <- mean(row, na.rm = TRUE)
    row
  }))
  
  # Title with spacing
  title_padding <- paste(rep("\n", title_margin_lines), collapse = "")
  plot_title    <- paste0("CV — Top ", top_n,
                          " most uncertain TEs in ", s, title_padding)
  
  annotation_columns <- data.frame(
    " "         = as.character(colData(se)$condition),
    row.names   = colnames(se),
    check.names = FALSE
  )
  
  outfile <- file.path(OUT_DIR, "plots",
                       paste0("04_heatmap_cv_top40_", s, ".pdf"))
  
  saved <- tryCatch({
    pdf(outfile, width = max(8, n_samples * 0.7 + 3), height = 10)
    pheatmap(
      mat_sub,
      color             = viridis(100),
      cluster_rows      = FALSE,
      cluster_cols      = FALSE,
      annotation_col    = annotation_columns,
      annotation_colors = condition_colors,
      show_rownames     = TRUE,
      show_colnames     = TRUE,
      fontsize_row      = 7,
      border_color      = cell_border_color,
      main              = plot_title,
      na_col            = "grey80",
      fontsize          = 10,
      main_fontface     = title_fontface,
      main_fontsize     = title_fontsize,
      silent            = FALSE
    )
    dev.off()
    cat("  ->", basename(outfile), "saved\n")
    TRUE
  }, error = function(e) {
    if (dev.cur() > 1) dev.off()
    cat("  ERROR for", s, ":", conditionMessage(e), "\n")
    FALSE
  })
}


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
  cat("  -> Plot 5 skipped: only one sample available.\n\n")
}


# ==============================================================================
# STEP 12: Final summary
# ==============================================================================

n_complete_te <- sum(complete.cases(mat_cv))

cat("=======================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=======================================================================\n")
cat("TEs analysed:              ", n_te,        "\n")
cat("Samples:                   ", n_samples,   "\n")
cat("Gibbs replicates:          ", n_gibbs,     "\n")
cat("Families found:            ", length(unique(copy_metrics$family)), "\n")
cat("TEs with complete CV data: ", n_complete_te,
    "(", round(n_complete_te / n_te * 100, 1), "% of total)\n")
cat("TEs using CV metric:       ", n_cv_regime,   "\n")
cat("TEs using Fano metric:     ", n_fano_regime, "\n")
cat("\nGlobal mean CV:            ",
    round(mean(copy_metrics$mean_cv, na.rm = TRUE), 4), "\n")
cat("Global median CV:          ",
    round(median(copy_metrics$mean_cv, na.rm = TRUE), 4), "\n")
cat("Most uncertain TE:         ", copy_metrics$te_copy[1],
    "(CV =", round(copy_metrics$mean_cv[1], 4), ")\n")
cat("\nOutput files saved in:", OUT_DIR, "\n")
cat("  tables/\n")
cat("    uncertainty_by_copy.tsv\n")
cat("    uncertainty_by_family.tsv\n")
cat("  plots/\n")
cat("    01_uncertainty_distribution.pdf\n")
cat("    02_abundance_vs_uncertainty.pdf\n")
cat("    03_cv_by_family_boxplot.pdf\n")
cat("    04_heatmap_cv_top40_<sample>.pdf  (one per sample)\n")
if (n_samples > 1) cat("    05_cv_violin_per_sample.pdf\n")
cat("=======================================================================\n")

sessionInfo()
