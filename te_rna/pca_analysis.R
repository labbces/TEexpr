# PCA ANALYSIS — DROUGHT STRESS TRANSCRIPTOMICS
# Groups: control | drought | rehydration (rehydration merged with drought)

# ============================================================================
# DEPENDENCIES:
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(patchwork)
library(DESeq2)
library(tximport)
library(ggrepel)

# ============================================================================
# COMMAND-LINE ARGUMENTS
# Cluster: Rscript pca_drought.R --species Sbicolor
# Interactive (RStudio): species_tag defined below is used as fallback
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else NULL
}

cli_species    <- parse_arg(args, "--species")
cli_salmon     <- parse_arg(args, "--salmon_root")
cli_metadata   <- parse_arg(args, "--metadata")
cli_output     <- parse_arg(args, "--output_dir")

# ============================================================================
# PARAMETERS — change these for each species/run
# ============================================================================

# --- Species tag (used to find sample folders by prefix) ---
# Overridden by --species argument when submitted to cluster
species_tag <- if (!is.null(cli_species)) cli_species else "Shybrid"  # overridden by --species

# --- Root directory containing salmon result folders for this species ---
salmon_root       <- if (!is.null(cli_salmon))   cli_salmon   else "~/pca_test/"  # overridden by --salmon_root

# --- Metadata table (columns: plant_code, genotype, run, bioproject, description, treatment, tissue) ---
metadata_file     <- if (!is.null(cli_metadata)) cli_metadata else "~/pca_test/meta_hybrid.tsv"  # overridden by --metadata

# --- Output directory ---
output_dir        <- if (!is.null(cli_output))   cli_output   else file.path("~/PCA_hybrid", species_tag)  # overridden by --output_dir

# --- Analysis mode ---
# "multi"  -> all experiments combined (colour = treatment, shape = genotype, label = plant_code)
# "single" -> one experiment only      (colour = treatment, shape = genotype)
# "both"   -> run both and save side-by-side panel
mode              <- "both"

# --- For single-experiment mode: which bioproject to use? ---
single_bioproject <- "PRJNA882367"

# --- CV filter threshold (keep genes with CV >= this value) ---
cv_threshold      <- 0.15

# --- Number of top variable genes fed into PCA ---
pca_top_genes     <- 1000

# --- Merge rehydration into drought? ---
merge_rehydration <- TRUE

# --- Colour palettes ---
# Multi mode: one hex colour per genotype, in the order they appear in the metadata.
# If NULL, falls back to viridis. Add or remove colours to match your number of genotypes.
colours_multi <- c("#bf7e46", "#ba4134", "#949941", "#ddc5a9",
                   "#5d556a", "#592E29", "#b03060", "#1a657a")

# Single mode: named vector with exactly "Control" and "Drought"
colours_single <- c(Control = "#f5c74d", Drought = "#8894b7")

# --- Títulos dos plots ---
title_multi  <- sprintf("%s — Múltiplos experimentos\n(cor = tratamento, forma = genótipo, rótulo = experimento)", species_tag)
title_single <- sprintf("%s — %s\n(cor = tratamento, forma = genótipo)", species_tag, single_bioproject)

# --- Nomes das legendas ---
legend_treatment  <- "Tratamento"
legend_genotype   <- "Genótipo"
legend_experiment <- "Experimento"

# --- Rótulos de tratamento (aparecem na legenda) ---
labels_treatment <- c(Control = "Controle", Drought = "Seca")

# ============================================================================
# HELPERS
# ============================================================================

# Build a hybrid tx2gene from quant.sf IDs
build_tx2gene <- function(sf_files) {
  
  cat("  [tx2gene] Reading feature IDs from quant.sf...\n")
  sf_names <- read.table(sf_files[1], header = TRUE, sep = "\t")$Name
  
  # -- 1. Separate genes and TEs -----------------------------------------------
  is_te    <- startsWith(sf_names, "TE_")
  gene_ids <- sf_names[!is_te]
  te_ids   <- sf_names[is_te]
  
  # -- 2. GENES: remove the isoform suffix (.1, .2, ...) ----------------------
  gene_rows <- data.frame(
    tx_id   = gene_ids,
    gene_id = sub("\\.\\d+$", "", gene_ids),
    stringsAsFactors = FALSE
  ) %>% distinct()
  
  cat(sprintf("  [tx2gene] Genes: %d transcripts -> %d unique genes\n",
              nrow(gene_rows), length(unique(gene_rows$gene_id))))
  
  # -- 3. TEs: remove "|coords|strand" then "_copyNNNN" -----------------------
  if (length(te_ids) == 0) {
    cat("  [tx2gene] No TE IDs found -- skipping TE mapping\n")
    te_rows <- data.frame(tx_id = character(0), gene_id = character(0),
                          stringsAsFactors = FALSE)
  } else {
    te_base <- sub("\\|.*", "", te_ids)
    te_repr <- sub("_copy.*", "", te_base)
    
    te_rows <- data.frame(tx_id   = te_ids,
                          gene_id = te_repr,
                          stringsAsFactors = FALSE) %>%
      distinct()
    
    cat(sprintf("  [tx2gene] TEs: %d instances -> %d representative elements\n",
                nrow(te_rows), length(unique(te_rows$gene_id))))
  }
  
  # -- 4. Combine and return ---------------------------------------------------
  tx2gene <- bind_rows(gene_rows, te_rows)
  cat(sprintf("  [tx2gene] Total: %d entries (%d gene transcripts + %d TE instances)\n",
              nrow(tx2gene), nrow(gene_rows), nrow(te_rows)))
  tx2gene
}


# Find all salmon quant.sf files for a species under salmon_root
find_sf_files <- function(salmon_root, species_tag) {
  species_dirs <- list.dirs(salmon_root, recursive = FALSE, full.names = TRUE)
  species_dirs <- species_dirs[startsWith(basename(species_dirs), species_tag)]
  
  if (length(species_dirs) == 0)
    stop(sprintf("No directories starting with '%s' found under:\n  %s", species_tag, salmon_root))
  
  cat(sprintf("  Found %d species director%s matching '%s'\n",
              length(species_dirs), ifelse(length(species_dirs) == 1, "y", "ies"), species_tag))
  
  all_sf <- c()
  for (sp_dir in species_dirs) {
    sample_dirs <- list.dirs(sp_dir, recursive = FALSE, full.names = TRUE)
    sf_paths    <- file.path(sample_dirs, "quant.sf")
    found       <- sf_paths[file.exists(sf_paths)]
    all_sf      <- c(all_sf, found)
  }
  
  if (length(all_sf) == 0)
    stop("No quant.sf files found. Check salmon_root and species_tag.")
  
  names(all_sf) <- basename(dirname(all_sf))
  cat(sprintf("  Found %d quant.sf files\n", length(all_sf)))
  all_sf
}


# Load and prepare metadata; merge rehydration -> drought if requested
load_metadata <- function(metadata_file, species_tag, merge_rehydration) {
  meta <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  required_cols <- c("plant_code", "genotype", "run", "bioproject", "treatment")
  missing_cols  <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0)
    stop(sprintf("Metadata missing columns: %s", paste(missing_cols, collapse = ", ")))
  
  meta$treatment <- tolower(trimws(meta$treatment))
  
  if (merge_rehydration) {
    n_rehy <- sum(meta$treatment %in% c("rehydration", "rehydrated", "rehidratacao", "rehidratación"))
    if (n_rehy > 0) {
      cat(sprintf("  Merging %d rehydration samples -> 'drought'\n", n_rehy))
      meta$treatment[meta$treatment %in% c("rehydration", "rehydrated",
                                           "rehidratacao", "rehidratación")] <- "drought"
    }
  }
  
  meta$treatment <- factor(meta$treatment,
                           levels = c("control", "drought"),
                           labels = c("Control", "Drought"))
  meta
}


# Coefficient-of-variation filter on raw count matrix
filter_by_cv <- function(counts_mat, cv_threshold) {
  gene_means <- rowMeans(counts_mat)
  gene_sds   <- apply(counts_mat, 1, sd)
  cv         <- ifelse(gene_means == 0, 0, gene_sds / gene_means)
  keep       <- cv >= cv_threshold
  cat(sprintf("  CV filter (>= %.0f%%): keeping %d / %d genes\n",
              cv_threshold * 100, sum(keep), length(keep)))
  counts_mat[keep, , drop = FALSE]
}


# Core PCA builder — returns a ggplot object
# @param txi       tximport object
# @param coldata   data.frame with rownames = run IDs
# @param color_col column name for point colour
# @param shape_col column name for point shape (NULL = no shape aesthetic)
# @param title     plot title
# @param colours   colour vector to use (passed from parameters)
build_pca_plot <- function(txi, coldata, color_col, shape_col = NULL,
                           title = "", colours = NULL, label_col = NULL) {
  
  design_formula <- as.formula(paste("~", color_col))
  
  dds <- DESeqDataSetFromTximport(
    txi     = txi,
    colData = coldata,
    design  = design_formula
  )
  
  raw_filtered <- filter_by_cv(counts(dds), cv_threshold)
  dds          <- dds[rownames(raw_filtered), ]
  dds          <- estimateSizeFactors(dds)
  
  vst      <- varianceStabilizingTransformation(dds, blind = TRUE)
  pca_data <- plotPCA(vst, intgroup = c(color_col, shape_col),
                      returnData = TRUE, ntop = pca_top_genes)
  pct_var  <- round(100 * attr(pca_data, "percentVar"))
  
  # Filled shape values (15=square, 16=circle, 17=triangle, 18=diamond, ...)
  filled_shapes <- c(16, 15, 17, 18, 21, 22, 23, 24, 25)
  
  if (!is.null(shape_col)) {
    pca_data[[shape_col]] <- as.factor(pca_data[[shape_col]])
    n_shapes <- nlevels(pca_data[[shape_col]])
    if (n_shapes > length(filled_shapes))
      warning(sprintf("shape_col '%s' has %d levels — consider using colour only.", shape_col, n_shapes))
    shape_values <- filled_shapes[seq_len(min(n_shapes, length(filled_shapes)))]
  }
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                            colour = .data[[color_col]],
                            shape  = if (!is.null(shape_col)) .data[[shape_col]] else NULL)) +
    geom_point(size = 4, alpha = 0.85, stroke = 0.5) +
    { if (!is.null(label_col))
      ggrepel::geom_text_repel(aes(label = .data[[label_col]]),
                               size = 3, max.overlaps = 20,
                               box.padding = 0.4, show.legend = FALSE)
      else NULL } +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle(title) +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 16, hjust = 0.5,
                                      margin = margin(b = 8)),
      legend.title     = element_text(face = "bold", size = 12),
      legend.text      = element_text(size = 11),
      legend.position  = "right",
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.title       = element_text(size = 12),
      axis.text        = element_text(size = 11),
      plot.background  = element_rect(fill = "white", colour = NA)
    )
  
  # Colour scale
  legend_name <- switch(color_col,
                        "treatment" = legend_treatment,
                        "genotype"  = legend_genotype,
                        color_col
  )
  if (color_col == "treatment") {
    p <- p + scale_colour_manual(values = colours, labels = labels_treatment, name = legend_name)
  } else if (!is.null(colours)) {
    p <- p + scale_colour_manual(values = colours, name = legend_name)
  } else {
    p <- p + scale_colour_viridis_d(option = "D", name = legend_name)
  }
  
  # Shape scale
  if (!is.null(shape_col)) {
    shape_name <- switch(shape_col,
                         "genotype"   = legend_genotype,
                         "plant_code" = legend_experiment,
                         shape_col
    )
    p <- p + scale_shape_manual(values = shape_values, name = shape_name)
  }
  
  p
}

# ============================================================================
# MAIN
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat(sprintf("PCA ANALYSIS — Species: %s\n", species_tag))
cat(sprintf("Mode: %s\n", mode))
cat("================================================================\n")

# --- 1. Find quant.sf files ---
cat("\n[1] Locating quant.sf files...\n")
all_sf <- find_sf_files(salmon_root, species_tag)

# --- 2. Load metadata ---
cat("\n[2] Loading metadata...\n")
meta <- load_metadata(metadata_file, species_tag, merge_rehydration)

common_runs <- intersect(names(all_sf), meta$run)
if (length(common_runs) == 0)
  stop("No overlap between quant.sf folder names and metadata 'run' column. Check run IDs.")

cat(sprintf("  Matched %d samples (metadata: %d | sf files: %d)\n",
            length(common_runs), nrow(meta), length(all_sf)))

sf_matched   <- all_sf[common_runs]
meta_matched <- meta %>%
  filter(run %in% common_runs) %>%
  arrange(match(run, common_runs))
rownames(meta_matched) <- meta_matched$run

na_treat <- is.na(meta_matched$treatment)
if (any(na_treat)) {
  cat(sprintf("  Warning: Dropping %d samples with unrecognised treatment\n", sum(na_treat)))
  meta_matched <- meta_matched[!na_treat, ]
  sf_matched   <- sf_matched[rownames(meta_matched)]
}

# --- 3. Build tx2gene (genes + TEs, entirely from .sf IDs) ---
cat("\n[3] Building tx2gene...\n")
tx2gene <- build_tx2gene(sf_matched)

# --- 4. tximport (shared for both modes) ---
cat("\n[4] Importing Salmon quantifications...\n")
txi_all <- tximport(
  files           = sf_matched,
  type            = "salmon",
  tx2gene         = tx2gene,
  ignoreTxVersion = TRUE,
  ignoreAfterBar  = TRUE
)
cat(sprintf("  Imported | Samples: %d | Genes (pre-filter): %d\n",
            ncol(txi_all$counts), nrow(txi_all$counts)))

plots <- list()

# ============================================================
# MODE: MULTI-EXPERIMENT
# colour = genotype, shape = plant_code (experiment code)
# ============================================================
if (mode %in% c("multi", "both")) {
  cat("\n[MULTI] Building multi-experiment PCA...\n")
  
  coldata_multi            <- meta_matched
  coldata_multi$treatment  <- factor(coldata_multi$treatment,
                                     levels = c("Control", "Drought"))
  coldata_multi$genotype   <- as.factor(coldata_multi$genotype)
  coldata_multi$plant_code <- as.factor(coldata_multi$plant_code)
  
  plots[["multi"]] <- build_pca_plot(
    txi       = txi_all,
    coldata   = coldata_multi,
    color_col = "treatment",
    shape_col = "genotype",
    label_col = "plant_code",
    title     = title_multi,
    colours   = colours_single
  )
  cat("  Multi-experiment PCA done\n")
}

# ============================================================
# MODE: SINGLE-EXPERIMENT
# colour = treatment, shape = genotype
# ============================================================
if (mode %in% c("single", "both")) {
  cat("\n[SINGLE] Building single-experiment PCA...\n")
  
  available_bp <- unique(meta_matched$bioproject)
  
  if (is.null(single_bioproject)) {
    cat(sprintf("  Available bioprojects: %s\n", paste(available_bp, collapse = ", ")))
    stop("Set 'single_bioproject' in the PARAMETERS section to one of the values above.")
  }
  
  if (!single_bioproject %in% available_bp)
    stop(sprintf("single_bioproject '%s' not found. Available: %s",
                 single_bioproject, paste(available_bp, collapse = ", ")))
  
  meta_single <- meta_matched %>% filter(bioproject == single_bioproject)
  cat(sprintf("  Bioproject: %s | Samples: %d\n", single_bioproject, nrow(meta_single)))
  
  idx        <- which(colnames(txi_all$counts) %in% meta_single$run)
  txi_single <- list(
    counts              = txi_all$counts[, idx, drop = FALSE],
    abundance           = txi_all$abundance[, idx, drop = FALSE],
    length              = txi_all$length[, idx, drop = FALSE],
    countsFromAbundance = txi_all$countsFromAbundance
  )
  class(txi_single) <- "list"
  
  coldata_single           <- meta_single
  coldata_single$treatment <- factor(coldata_single$treatment,
                                     levels = c("Control", "Drought"))
  coldata_single$genotype  <- as.factor(coldata_single$genotype)
  
  plots[["single"]] <- build_pca_plot(
    txi       = txi_single,
    coldata   = coldata_single,
    color_col = "treatment",
    shape_col = "genotype",
    title     = title_single,
    colours   = colours_single
  )
  cat("  Single-experiment PCA done\n")
}

# ============================================================
# SAVE OUTPUT
# ============================================================
cat("\n[5] Saving output...\n")

if (mode == "both") {
  panel <- plots[["multi"]] + plots[["single"]] +
    plot_annotation(
      tag_levels = "A",
      theme = theme(
        plot.tag        = element_text(size = 18, face = "bold"),
        plot.background = element_rect(fill = "white", colour = NA)
      )
    ) +
    plot_layout(ncol = 2, widths = c(1, 1))
  
  out_file <- file.path(output_dir, sprintf("pca_panel_%s.png", species_tag))
  ggsave(out_file, plot = panel, width = 40, height = 18, units = "cm", dpi = 320, bg = "white")
  cat(sprintf("  Panel saved: %s\n", out_file))
  
} else {
  out_file <- file.path(output_dir, sprintf("pca_%s_%s.png", mode, species_tag))
  ggsave(out_file, plot = plots[[mode]], width = 20, height = 18, units = "cm", dpi = 320, bg = "white")
  cat(sprintf("  Plot saved: %s\n", out_file))
}

cat("\n================================================================\n")
cat("DONE\n")
cat("================================================================\n")
