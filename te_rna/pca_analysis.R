# ============================================================================
# PCA ANALYSIS — DROUGHT STRESS TRANSCRIPTOMICS
# Supports: multi-experiment mode and single-experiment mode
# Groups: control | drought | rehydration (rehydration merged with drought)
# ============================================================================
# DEPENDENCIES
# ============================================================================

library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
library(DESeq2)
library(tximport)
library(ggrepel)

# ============================================================================
# COMMAND-LINE ARGUMENTS
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
cli_multi_runs <- parse_arg(args, "--multi_runs")        # path to .txt with one SRR per line
cli_single_bp  <- parse_arg(args, "--single_bioproject")

# ============================================================================
# PARAMETERS
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
# "multi"  -> all experiments combined (colour = genotype, shape = treatment, label = plant_code)
# "single" -> one experiment only      (colour = genotype, shape = treatment)
# "both"   -> run both and save side-by-side panel
mode              <- "both"

# --- For single-experiment mode: which bioproject to use? ---
single_bioproject <- if (!is.null(cli_single_bp)) cli_single_bp else "PRJNA882367"  # overridden by --single_bioproject

# --- For multi-experiment mode: optional whitelist of SRRs (one per line in a .txt) ---
# If provided, ONLY these runs are used in the multi PCA. If NULL, all available runs are used.
multi_runs_file   <- if (!is.null(cli_multi_runs)) cli_multi_runs else NULL  # set via --multi_runs

# --- CV filter threshold (keep genes with CV >= this value) ---
cv_threshold      <- 0.15

# --- Number of top variable genes fed into PCA ---
pca_top_genes     <- 1000

# --- Merge rehydration into drought? ---
merge_rehydration <- TRUE

# --- Colour palette (msc_palette — 40 colours) ---
# Organised in 10 hue families. If NULL, falls back to viridis automatically.
colours_multi <- c(
  # Dark Mahogany
  "#300B07", "#481B17", "#592E29", "#734C48",
  # Crimson Clay
  "#6A0A00", "#952418", "#ba4134", "#D8675B",
  # Caramel
  "#4C2300", "#723908", "#9E5D25", "#bf7e46",
  # Saffron
  "#936805", "#B4851A", "#DBA938", "#f5c75d",
  # Olive Drab
  "#52560B", "#787D24", "#949941", "#B6BB64",
  # Verde/Ciano
  "#5C8C3A", "#3A8C6A", "#3A8C8C", "#3A6A8C",
  # Teal Blue
  "#023A4A", "#0D5063", "#1a657a", "#317689",
  # Steel Blue
  "#0A1842", "#38497C", "#5E6D9A", "#8894b7",
  # Dusty Plum
  "#1F1631", "#2A2238", "#433C50", "#5d556a",
  # Raspberry
  "#2D0011", "#4F001D", "#70032B", "#b03060"
)

# --- Anchor colours — one per family, used first for small genotype counts ---
# With <= 10 genotypes: only anchors are used (maximum contrast).
# With > 10 genotypes: full palette is sampled evenly.
colours_anchors <- c(
  "#592E29",  # Dark Mahogany
  "#ba4134",  # Crimson Clay
  "#bf7e46",  # Caramel
  "#f5c75d",  # Saffron
  "#949941",  # Olive Drab
  "#5C8C3A",  # Verde/Ciano
  "#1a657a",  # Teal Blue
  "#8894b7",  # Steel Blue
  "#5d556a",  # Dusty Plum
  "#b03060"   # Raspberry
)

# --- Treatment shapes ---
# 16 = circle (filled), 17 = triangle (filled)
shapes_treatment <- c(Control = 16, Drought = 17)

# --- Plot titles ---
title_multi  <- sprintf("%s — Multi-experiment", species_tag)
# title_single is computed dynamically in the SINGLE block below

# --- Nomes das legendas ---
legend_treatment  <- "Treatment"
legend_genotype   <- "Genotype"
legend_experiment <- "Experiment"

# --- Treatment labels (shown in legend) ---
labels_treatment <- c(Control = "Control", Drought = "Drought")

# ============================================================================
# HELPERS
# ============================================================================

#' Select colours for n genotypes:
#' - if n <= 10: sample evenly from anchors (one per family, max contrast)
#' - if n > 10:  sample evenly from full palette
select_colours <- function(n, anchors = colours_anchors, palette = colours_multi) {
  if (is.null(palette) && is.null(anchors)) return(NULL)
  if (n <= length(anchors)) {
    idx <- round(seq(1, length(anchors), length.out = n))
    return(anchors[idx])
  } else if (!is.null(palette) && length(palette) >= n) {
    idx <- round(seq(1, length(palette), length.out = n))
    return(palette[idx])
  } else {
    warning(sprintf("Not enough colours (%d) for %d genotypes — falling back to viridis.",
                    length(palette), n))
    return(NULL)
  }
}

#' Build a hybrid tx2gene entirely from quant.sf IDs — no GTF needed.
#'
#' Gene transcripts : e.g. Sobic.010G000100.1  ->  remove trailing ".N"  ->  Sobic.010G000100
#' TE instances     : e.g. TE_00007302_copy0001|NC_008360.1:400945-401016|-
#'                        ->  remove "|..." then "_copy..."  ->  TE_00007302
#'
#' Detection rules: starts with "TE_" -> TE; otherwise -> gene transcript
#'
#' @param sf_files  named character vector of quant.sf paths (from find_sf_files)
#' @return data.frame with columns tx_id, gene_id  (ready for tximport)
build_tx2gene <- function(sf_files) {

  cat("  [tx2gene] Reading feature IDs from quant.sf...\n")
  sf_names <- read.table(sf_files[1], header = TRUE, sep = "\t")$Name

  # -- 1. Separate genes and TEs -----------------------------------------------
  is_te    <- startsWith(sf_names, "TE_")
  gene_ids <- sf_names[!is_te]
  te_ids   <- sf_names[is_te]

  # -- 2. GENES: remove the isoform suffix (.1, .2, ...) ----------------------
  # Remove isoform suffix — handles multiple patterns:
  #   Sobic.010G000100.1  -> Sobic.010G000100  (suffix: .N)
  #   Sof_g2.t1           -> Sof_g2            (suffix: .tN)
  #   Zm00001eb000010_T001 -> Zm00001eb000010   (suffix: _TN or _tN)
  gene_id_clean <- gene_ids
  gene_id_clean <- sub("\\.[Tt]\\d+$", "", gene_id_clean)   # remove .t1, .T1
  gene_id_clean <- sub("\\.\\d+$",     "", gene_id_clean)   # remove .1, .2
  gene_id_clean <- sub("_[Tt]\\d+$",    "", gene_id_clean)   # remove _T001, _t1

  gene_rows <- data.frame(
    tx_id   = gene_ids,
    gene_id = gene_id_clean,
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


#' Find all salmon quant.sf files for a species under salmon_root
find_sf_files <- function(salmon_root, species_tag) {
  species_dirs <- list.dirs(salmon_root, recursive = FALSE, full.names = TRUE)
  species_dirs <- species_dirs[startsWith(basename(species_dirs), species_tag)]

  if (length(species_dirs) == 0)
    stop(sprintf("No directories starting with '%s' found under:\n  %s", species_tag, salmon_root))

  cat(sprintf("  Found %d species director%s matching '%s'\n",
              length(species_dirs), ifelse(length(species_dirs) == 1, "y", "ies"), species_tag))

  all_sf <- c()
  for (sp_dir in species_dirs) {
    # Structure: <sp_dir>/salmon_results/<SRR>/quant.sf
    # Support both structures:
    #   cluster: <sp_dir>/salmon_results/<SRR>/quant.sf
    #   local:   <sp_dir>/<SRR>/quant.sf
    results_dir <- file.path(sp_dir, "salmon_results")
    if (dir.exists(results_dir)) {
      sample_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
    } else {
      cat(sprintf("  Note: no salmon_results subdir, reading directly from %s\n", sp_dir))
      sample_dirs <- list.dirs(sp_dir, recursive = FALSE, full.names = TRUE)
    }
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


#' Load and prepare metadata; merge rehydration -> drought if requested
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


#' Coefficient-of-variation filter on raw count matrix
#' Reports genes vs TEs separately in log output
filter_by_cv <- function(counts_mat, cv_threshold) {
  gene_means <- rowMeans(counts_mat)
  gene_sds   <- apply(counts_mat, 1, sd)
  cv         <- ifelse(gene_means == 0, 0, gene_sds / gene_means)
  keep       <- cv >= cv_threshold

  # Separate genes and TEs by ID prefix
  all_ids    <- rownames(counts_mat)
  is_te      <- startsWith(all_ids, "TE_")
  kept_ids   <- all_ids[keep]
  kept_te    <- startsWith(kept_ids, "TE_")

  n_genes_in  <- sum(!is_te);  n_genes_out <- sum(!kept_te)
  n_te_in     <- sum(is_te);   n_te_out    <- sum(kept_te)

  cat(sprintf("  CV filter (>= %.0f%%):\n", cv_threshold * 100))
  cat(sprintf("    Genes: %d / %d kept (%.1f%%)\n",
              n_genes_out, n_genes_in, 100 * n_genes_out / max(n_genes_in, 1)))
  cat(sprintf("    TEs:   %d / %d kept (%.1f%%)\n",
              n_te_out,    n_te_in,    100 * n_te_out    / max(n_te_in, 1)))
  cat(sprintf("    Total: %d / %d kept\n", sum(keep), length(keep)))

  counts_mat[keep, , drop = FALSE]
}


#' Core PCA builder — returns a ggplot object
#' @param txi       tximport object
#' @param coldata   data.frame with rownames = run IDs
#' @param color_col column name for point colour
#' @param shape_col column name for point shape (NULL = no shape aesthetic)
#' @param title     plot title
#' @param colours   colour vector to use (passed from parameters)
build_pca_plot <- function(txi, coldata, color_col, shape_col = NULL,
                           title = "", colours = NULL, label_col = NULL,
                           subtitle = NULL) {

  design_formula <- as.formula(paste("~", color_col))

  dds <- DESeqDataSetFromTximport(
    txi     = txi,
    colData = coldata,
    design  = design_formula
  )

  raw_filtered <- filter_by_cv(counts(dds), cv_threshold)
  dds          <- dds[rownames(raw_filtered), ]
  dds          <- estimateSizeFactors(dds)

  # Build gene/TE count subtitle
  kept_ids      <- rownames(raw_filtered)
  n_genes_kept  <- sum(!startsWith(kept_ids, "TE_"))
  n_te_kept     <- sum(startsWith(kept_ids, "TE_"))
  cv_subtitle   <- sprintf("Genes: %d | TEs: %d (CV >= %.0f%%)",
                            n_genes_kept, n_te_kept, cv_threshold * 100)
  if (is.null(subtitle)) subtitle <- cv_subtitle

  vst         <- varianceStabilizingTransformation(dds, blind = TRUE)
  intgroup    <- unique(c(color_col, shape_col, label_col))  # include label_col if provided
  intgroup    <- intgroup[!is.na(intgroup) & intgroup %in% colnames(colData(vst))]
  pca_data    <- plotPCA(vst, intgroup = intgroup,
                         returnData = TRUE, ntop = pca_top_genes)
  pct_var  <- round(100 * attr(pca_data, "percentVar"))

  if (!is.null(shape_col)) {
    pca_data[[shape_col]] <- as.factor(pca_data[[shape_col]])
    n_shapes <- nlevels(pca_data[[shape_col]])
    if (shape_col == "treatment") {
      # Fixed shapes: circle = Control, triangle = Drought
      shape_values <- shapes_treatment
    } else {
      filled_shapes <- c(16, 15, 17, 18, 21, 22, 23, 24, 25)
      if (n_shapes > length(filled_shapes))
        warning(sprintf("shape_col '%s' has %d levels — consider using colour only.", shape_col, n_shapes))
      shape_values <- filled_shapes[seq_len(min(n_shapes, length(filled_shapes)))]
    }
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
    { if (!is.null(subtitle)) labs(subtitle = subtitle) else NULL } +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 16, hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = 11, hjust = 0.5, colour = "grey40",
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
  if (!is.null(colours)) {
    p <- p + scale_colour_manual(values = colours, name = legend_name)
  } else {
    p <- p + scale_colour_viridis_d(option = "D", name = legend_name)
  }

  # Shape scale
  if (!is.null(shape_col)) {
    shape_name <- switch(shape_col,
      "treatment"  = legend_treatment,
      "genotype"   = legend_genotype,
      "plant_code" = legend_experiment,
      shape_col
    )
    if (shape_col == "treatment") {
      p <- p + scale_shape_manual(values = shape_values, labels = labels_treatment, name = shape_name)
    } else {
      p <- p + scale_shape_manual(values = shape_values, name = shape_name)
    }
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
  ignoreTxVersion = FALSE,
  ignoreAfterBar  = FALSE
)
cat(sprintf("  Imported | Samples: %d | Genes (pre-filter): %d\n",
            ncol(txi_all$counts), nrow(txi_all$counts)))

plots <- list()

# ============================================================
# MODE: MULTI-EXPERIMENT
# colour = genotype, shape = treatment, label = plant_code
# ============================================================
if (mode %in% c("multi", "both")) {
  cat("\n[MULTI] Building multi-experiment PCA...\n")

  # Optional whitelist of SRRs for the multi analysis
  if (!is.null(multi_runs_file)) {
    if (!file.exists(multi_runs_file))
      stop(sprintf("multi_runs file not found: %s", multi_runs_file))

    multi_runs <- readLines(multi_runs_file)
    multi_runs <- trimws(multi_runs)
    multi_runs <- multi_runs[nchar(multi_runs) > 0]   # drop blank lines

    cat(sprintf("  Whitelist provided (%s): %d SRRs requested\n",
                multi_runs_file, length(multi_runs)))

    runs_found   <- intersect(multi_runs, rownames(meta_matched))
    runs_missing <- setdiff(multi_runs, rownames(meta_matched))

    if (length(runs_missing) > 0)
      cat(sprintf("  Warning: %d whitelisted SRRs not found in data and skipped: %s\n",
                  length(runs_missing), paste(runs_missing, collapse = ", ")))

    if (length(runs_found) == 0)
      stop("None of the whitelisted SRRs were found in the available data.")

    cat(sprintf("  Using %d SRRs for multi-experiment PCA\n", length(runs_found)))

    meta_m <- meta_matched[meta_matched$run %in% runs_found, , drop = FALSE]
    meta_m <- meta_m[order(match(meta_m$run, colnames(txi_all$counts))), ]
    idx_m  <- match(meta_m$run, colnames(txi_all$counts))
    txi_m  <- list(
      counts              = txi_all$counts[, idx_m, drop = FALSE],
      abundance           = txi_all$abundance[, idx_m, drop = FALSE],
      length              = txi_all$length[, idx_m, drop = FALSE],
      countsFromAbundance = txi_all$countsFromAbundance
    )
    class(txi_m) <- "list"
  } else {
    cat(sprintf("  No whitelist provided — using all %d available SRRs\n", nrow(meta_matched)))
    meta_m <- meta_matched[order(match(meta_matched$run, colnames(txi_all$counts))), ]
    txi_m  <- txi_all
  }

  coldata_multi            <- meta_m
  rownames(coldata_multi)  <- coldata_multi$run
  coldata_multi$treatment  <- factor(coldata_multi$treatment,
                                      levels = c("Control", "Drought"))
  coldata_multi$genotype   <- as.factor(coldata_multi$genotype)
  coldata_multi$plant_code <- as.factor(coldata_multi$plant_code)

  # Build named colour vector — anchors first, full palette if > 10 genotypes
  geno_levels         <- levels(coldata_multi$genotype)
  n_geno              <- length(geno_levels)
  selected            <- select_colours(n_geno)
  colours_multi_named <- if (!is.null(selected)) setNames(selected, geno_levels) else NULL

  plots[["multi"]] <- build_pca_plot(
    txi       = txi_m,
    coldata   = coldata_multi,
    color_col = "genotype",
    shape_col = "treatment",
    label_col = "plant_code",
    title     = title_multi,
    colours   = colours_multi_named
  )
  cat("  Multi-experiment PCA done\n")
}

# ============================================================
# MODE: SINGLE-EXPERIMENT
# colour = genotype, shape = treatment
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

  # Use plant_codes from this experiment for the title (matching panel A labels)
  plant_codes_single <- paste(sort(unique(meta_single$plant_code)), collapse = ", ")
  title_single       <- sprintf("%s — %s", species_tag, plant_codes_single)

  # Align coldata and txi to same samples in same order
  runs_single <- intersect(meta_single$run, colnames(txi_all$counts))
  meta_single <- meta_single[meta_single$run %in% runs_single, , drop = FALSE]
  meta_single <- meta_single[order(match(meta_single$run, colnames(txi_all$counts))), ]

  idx        <- match(meta_single$run, colnames(txi_all$counts))
  txi_single <- list(
    counts              = txi_all$counts[, idx, drop = FALSE],
    abundance           = txi_all$abundance[, idx, drop = FALSE],
    length              = txi_all$length[, idx, drop = FALSE],
    countsFromAbundance = txi_all$countsFromAbundance
  )
  class(txi_single) <- "list"

  coldata_single           <- meta_single
  rownames(coldata_single) <- coldata_single$run
  coldata_single$treatment <- factor(coldata_single$treatment,
                                      levels = c("Control", "Drought"))
  coldata_single$genotype  <- as.factor(coldata_single$genotype)

  # Build named colour vector for single-experiment genotypes
  geno_levels_s        <- levels(coldata_single$genotype)
  n_geno_s             <- length(geno_levels_s)
  selected_s           <- select_colours(n_geno_s)
  colours_single_named <- if (!is.null(selected_s)) setNames(selected_s, geno_levels_s) else NULL

  plots[["single"]] <- build_pca_plot(
    txi       = txi_single,
    coldata   = coldata_single,
    color_col = "genotype",
    shape_col = "treatment",
    title     = title_single,
    colours   = colours_single_named
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

# ============================================================
# SUMMARY TABLES
# ============================================================
cat("\n[6] Generating summary tables...\n")

# Table 1: filtering effect (genes and TEs before/after CV filter)
all_ids_total  <- rownames(txi_all$counts)
n_genes_total  <- sum(!startsWith(all_ids_total, "TE_"))
n_te_total     <- sum(startsWith(all_ids_total, "TE_"))

# Re-run filter to get post-filter counts (using meta_matched samples)
dds_summary    <- DESeqDataSetFromTximport(
  txi     = txi_all,
  colData = meta_matched[, c("treatment", "genotype", "plant_code"), drop = FALSE],
  design  = ~ treatment
)
dds_summary    <- estimateSizeFactors(dds_summary)
filt_summary   <- filter_by_cv(counts(dds_summary), cv_threshold)
kept_ids_summ  <- rownames(filt_summary)
n_genes_kept   <- sum(!startsWith(kept_ids_summ, "TE_"))
n_te_kept      <- sum(startsWith(kept_ids_summ, "TE_"))

filter_table <- data.frame(
  Feature   = c("Protein-coding genes", "Transposable elements (TEs)", "Total"),
  Before_CV = c(n_genes_total, n_te_total, n_genes_total + n_te_total),
  After_CV  = c(n_genes_kept,  n_te_kept,  n_genes_kept  + n_te_kept),
  Retained_pct = round(100 * c(n_genes_kept, n_te_kept, n_genes_kept + n_te_kept) /
                           c(n_genes_total, n_te_total, n_genes_total + n_te_total), 1)
)
colnames(filter_table)[4] <- "Retained (%)"

out_filter <- file.path(output_dir, sprintf("summary_cv_filter_%s.tsv", species_tag))
write.table(filter_table, file = out_filter, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Filtering summary saved: %s\n", out_filter))
cat("\n  CV Filtering Effect:\n")
print(filter_table, row.names = FALSE)

# Table 2: species summary (genotypes, experiments, tissues)
species_summary <- meta_matched %>%
  summarise(
    Species         = species_tag,
    N_samples       = n(),
    N_genotypes     = n_distinct(genotype),
    N_experiments   = n_distinct(plant_code),
    N_bioprojects   = n_distinct(bioproject),
    Tissues         = paste(sort(unique(tissue)), collapse = ", "),
    Treatments      = paste(sort(unique(as.character(treatment))), collapse = ", ")
  )

out_species <- file.path(output_dir, sprintf("summary_species_%s.tsv", species_tag))
write.table(species_summary, file = out_species, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  Species summary saved: %s\n", out_species))
cat("\n  Species Summary:\n")
print(as.data.frame(species_summary), row.names = FALSE)

cat("\n================================================================\n")
cat("DONE\n")
cat("================================================================\n")
