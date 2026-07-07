# ============================================================================
# RUV ANALYSIS — DROUGHT STRESS TRANSCRIPTOMICS
# Method: RUVr (residuals)
# Outputs: corrected counts + VST + W factors for each k
# ============================================================================

library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
library(ggrepel)
library(DESeq2)
library(tximport)
library(RUVSeq)
library(edgeR)

# ============================================================================
# COMMAND-LINE ARGUMENTS
# Cluster: Rscript ruv_analysis.R --species Sbicolor --salmon_root /path ...
# Interactive: fallback values in PARAMETERS section are used
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx < length(args)) args[idx + 1] else NULL
}

cli_species   <- parse_arg(args, "--species")
cli_salmon    <- parse_arg(args, "--salmon_root")
cli_metadata  <- parse_arg(args, "--metadata")
cli_output    <- parse_arg(args, "--output_dir")
cli_runs      <- parse_arg(args, "--runs_file")     # same whitelist as PCA multi
cli_single_bp <- parse_arg(args, "--single_bioproject")
cli_k_values  <- parse_arg(args, "--k_values")         # e.g. "6:11" or "1:5"

# ============================================================================
# PARAMETERS — change these for each species/run
# ============================================================================

species_tag   <- if (!is.null(cli_species))  cli_species  else "Shybrid"
salmon_root   <- if (!is.null(cli_salmon))   cli_salmon   else "~/pca_test/"
metadata_file <- if (!is.null(cli_metadata)) cli_metadata else "~/pca_test/meta_hybrid.tsv"
output_dir    <- if (!is.null(cli_output))   cli_output   else file.path("~/RUV", species_tag)
runs_file     <- if (!is.null(cli_runs))     cli_runs     else NULL  # path to SRR whitelist .txt

# --- RUV settings ---
# k_values: overridden via --k_values "6:11"
k_values      <- if (!is.null(cli_k_values)) eval(parse(text = cli_k_values)) else 1:5
cv_threshold      <- 0.5    # same CV filter as PCA
pca_top_genes     <- 1000   # top variable genes for PCA plots
n_negative_ctrl   <- 1000   # number of empirical negative controls for RUVg (lowest CV genes)

# --- Merge rehydration into drought? ---
merge_rehydration <- TRUE

# --- Plot aesthetics (same as PCA script) ---
shapes_treatment  <- c(Control = 16, Drought = 17)
legend_treatment  <- "Tratamento"
legend_genotype   <- "Genótipo"
legend_experiment <- "Experimento"
labels_treatment  <- c(Control = "Controle", Drought = "Seca")

# --- Colour palette per genotype (msc_palette — 40 colours) ---
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

# Helper: selects colours prioritising anchors for small n
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

# ============================================================================
# HELPERS (shared with PCA script)
# ============================================================================

build_tx2gene <- function(sf_files) {
  cat("  [tx2gene] Reading feature IDs from quant.sf...\n")
  sf_names <- read.table(sf_files[1], header = TRUE, sep = "\t")$Name
  
  is_te    <- startsWith(sf_names, "TE_")
  gene_ids <- sf_names[!is_te]
  te_ids   <- sf_names[is_te]
  
  # Remove isoform suffix — handles multiple patterns:
  #   Sobic.010G000100.1  -> Sobic.010G000100  (suffix: .N)
  #   Sof_g2.t1           -> Sof_g2            (suffix: .tN)
  #   Zm00001eb000010_T001 -> Zm00001eb000010  (suffix: _TN or _tN)
  gene_id_clean <- gene_ids
  gene_id_clean <- sub("\\.[Tt]\\d+$", "", gene_id_clean)
  gene_id_clean <- sub("\\.\\d+$",     "", gene_id_clean)
  gene_id_clean <- sub("_[Tt]\\d+$",  "", gene_id_clean)
  
  gene_rows <- data.frame(
    tx_id   = gene_ids,
    gene_id = gene_id_clean,
    stringsAsFactors = FALSE
  ) %>% distinct()
  
  cat(sprintf("  [tx2gene] Genes: %d transcripts -> %d unique genes\n",
              nrow(gene_rows), length(unique(gene_rows$gene_id))))
  
  if (length(te_ids) == 0) {
    cat("  [tx2gene] No TE IDs found -- skipping TE mapping\n")
    te_rows <- data.frame(tx_id = character(0), gene_id = character(0),
                          stringsAsFactors = FALSE)
  } else {
    te_base <- sub("\\|.*", "", te_ids)
    te_repr <- sub("_copy.*", "", te_base)
    te_rows <- data.frame(tx_id   = te_ids,
                          gene_id = te_repr,
                          stringsAsFactors = FALSE) %>% distinct()
    cat(sprintf("  [tx2gene] TEs: %d instances -> %d representative elements\n",
                nrow(te_rows), length(unique(te_rows$gene_id))))
  }
  
  tx2gene <- bind_rows(gene_rows, te_rows)
  cat(sprintf("  [tx2gene] Total: %d entries (%d gene transcripts + %d TE instances)\n",
              nrow(tx2gene), nrow(gene_rows), nrow(te_rows)))
  tx2gene
}


find_sf_files <- function(salmon_root, species_tag) {
  species_dirs <- list.dirs(salmon_root, recursive = FALSE, full.names = TRUE)
  species_dirs <- species_dirs[startsWith(basename(species_dirs), species_tag)]
  
  if (length(species_dirs) == 0)
    stop(sprintf("No directories starting with '%s' found under:\n  %s", species_tag, salmon_root))
  
  cat(sprintf("  Found %d species director%s matching '%s'\n",
              length(species_dirs), ifelse(length(species_dirs) == 1, "y", "ies"), species_tag))
  
  all_sf <- c()
  for (sp_dir in species_dirs) {
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


load_metadata <- function(metadata_file, merge_rehydration) {
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


filter_by_cv <- function(counts_mat, cv_threshold) {
  gene_means <- rowMeans(counts_mat)
  gene_sds   <- apply(counts_mat, 1, sd)
  cv         <- ifelse(gene_means == 0, 0, gene_sds / gene_means)
  keep       <- cv >= cv_threshold
  cat(sprintf("  CV filter (>= %.0f%%): keeping %d / %d genes\n",
              cv_threshold * 100, sum(keep), length(keep)))
  counts_mat[keep, , drop = FALSE]
}


#' Build PCA plot from a DESeq2 VST object
build_pca_plot <- function(vst_obj, coldata, title = "") {
  color_col <- "genotype"
  shape_col <- "treatment"
  label_col <- "plant_code"
  
  intgroup  <- unique(c(color_col, shape_col, label_col))
  intgroup  <- intgroup[!is.na(intgroup) & intgroup %in% colnames(colData(vst_obj))]
  pca_data  <- plotPCA(vst_obj, intgroup = intgroup,
                       returnData = TRUE, ntop = pca_top_genes)
  pct_var   <- round(100 * attr(pca_data, "percentVar"))
  
  # Build gene/TE count subtitle
  feat_ids     <- rownames(assay(vst_obj))
  n_genes_plot <- sum(!startsWith(feat_ids, "TE_"))
  n_te_plot    <- sum(startsWith(feat_ids, "TE_"))
  subtitle     <- sprintf("Genes: %d | TEs: %d (CV >= %.0f%%)",
                          n_genes_plot, n_te_plot, cv_threshold * 100)
  
  pca_data[[shape_col]] <- as.factor(pca_data[[shape_col]])
  pca_data[[color_col]] <- as.factor(pca_data[[color_col]])
  
  n_geno <- nlevels(pca_data[[color_col]])
  
  p <- ggplot(pca_data, aes(x = PC1, y = PC2,
                            colour = .data[[color_col]],
                            shape  = .data[[shape_col]])) +
    geom_point(size = 4, alpha = 0.85, stroke = 0.5) +
    ggrepel::geom_text_repel(aes(label = .data[[label_col]]),
                             size = 3, max.overlaps = 20,
                             box.padding = 0.4, show.legend = FALSE) +
    xlab(paste0("PC1: ", pct_var[1], "% variância")) +
    ylab(paste0("PC2: ", pct_var[2], "% variância")) +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    { n_col   <- nlevels(as.factor(pca_data[[color_col]]))
    sel_col <- select_colours(n_col)
    if (!is.null(sel_col))
      scale_colour_manual(
        values = setNames(sel_col, levels(as.factor(pca_data[[color_col]]))),
        name   = legend_genotype)
    else
      scale_colour_viridis_d(option = "D", name = legend_genotype) } +
    scale_shape_manual(values  = shapes_treatment,
                       labels  = labels_treatment,
                       name    = legend_treatment) +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5,
                                      margin = margin(b = 4)),
      plot.subtitle    = element_text(size = 10, hjust = 0.5, colour = "grey40",
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
  p
}

# ============================================================================
# MAIN
# ============================================================================

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("================================================================\n")
cat(sprintf("RUVr ANALYSIS — Species: %s\n", species_tag))
cat(sprintf("k values: %s\n", paste(k_values, collapse = ", ")))
cat("================================================================\n")

# --- 1. Find quant.sf files ---
cat("\n[1] Locating quant.sf files...\n")
all_sf <- find_sf_files(salmon_root, species_tag)

# --- 2. Load metadata ---
cat("\n[2] Loading metadata...\n")
meta <- load_metadata(metadata_file, merge_rehydration)

common_runs <- intersect(names(all_sf), meta$run)
if (length(common_runs) == 0)
  stop("No overlap between quant.sf folder names and metadata 'run' column.")

cat(sprintf("  Matched %d samples (metadata: %d | sf files: %d)\n",
            length(common_runs), nrow(meta), length(all_sf)))

sf_matched   <- all_sf[common_runs]
meta_matched <- meta %>%
  filter(run %in% common_runs) %>%
  arrange(match(run, common_runs))
rownames(meta_matched) <- meta_matched$run

# Drop samples with unrecognised treatment
na_treat <- is.na(meta_matched$treatment)
if (any(na_treat)) {
  cat(sprintf("  Warning: Dropping %d samples with unrecognised treatment\n", sum(na_treat)))
  meta_matched <- meta_matched[!na_treat, ]
  sf_matched   <- sf_matched[rownames(meta_matched)]
}

# --- 3. Apply SRR whitelist (leaf filter) ---
if (!is.null(runs_file)) {
  if (!file.exists(runs_file))
    stop(sprintf("runs_file not found: %s", runs_file))
  
  whitelist   <- trimws(readLines(runs_file))
  whitelist   <- whitelist[nchar(whitelist) > 0]
  runs_found  <- intersect(whitelist, rownames(meta_matched))
  runs_miss   <- setdiff(whitelist, rownames(meta_matched))
  
  cat(sprintf("  Whitelist: %d requested | %d found | %d missing\n",
              length(whitelist), length(runs_found), length(runs_miss)))
  if (length(runs_miss) > 0)
    cat(sprintf("  Missing: %s\n", paste(runs_miss, collapse = ", ")))
  if (length(runs_found) == 0)
    stop("None of the whitelisted SRRs found in data.")
  
  meta_matched <- meta_matched[runs_found, , drop = FALSE]
  meta_matched <- meta_matched[order(match(meta_matched$run, names(sf_matched))), ]
  sf_matched   <- sf_matched[meta_matched$run]
} else {
  cat(sprintf("  No whitelist — using all %d samples\n", nrow(meta_matched)))
}

# --- 4. Build tx2gene ---
cat("\n[3] Building tx2gene...\n")
tx2gene <- build_tx2gene(sf_matched)

# --- 5. tximport ---
cat("\n[4] Importing Salmon quantifications...\n")
txi <- tximport(
  files           = sf_matched,
  type            = "salmon",
  tx2gene         = tx2gene,
  ignoreTxVersion = FALSE,
  ignoreAfterBar  = FALSE
)
cat(sprintf("  Imported | Samples: %d | Genes (pre-filter): %d\n",
            ncol(txi$counts), nrow(txi$counts)))

# --- 6. Build DESeq2 object and filter ---
cat("\n[5] Building DESeq2 object...\n")

coldata           <- meta_matched
coldata$treatment <- factor(coldata$treatment, levels = c("Control", "Drought"))
coldata$genotype  <- as.factor(coldata$genotype)
coldata$plant_code <- as.factor(coldata$plant_code)

dds_raw <- DESeqDataSetFromTximport(
  txi     = txi,
  colData = coldata,
  design  = ~ treatment
)

# Filter: keep only expressed genes
keep <- rowSums(counts(dds_raw)) > 0
dds_raw <- dds_raw[keep, ]
head(as.data.frame(assays(dds_raw)))
cat(sprintf("  After expression filter: %d genes retained\n", sum(keep)))

# CV filter
dds_raw_pre  <- estimateSizeFactors(dds_raw) #remove sequencing effort effect
raw_filtered <- filter_by_cv(counts(dds_raw_pre), cv_threshold)
dds_raw      <- dds_raw[rownames(raw_filtered), ]

# --- 7. Raw PCA (before correction — control panel) ---
cat("\n[6] Building raw PCA (pre-correction)...\n")

dds_raw_vst <- dds_raw
dds_raw_vst <- estimateSizeFactors(dds_raw_vst)
vst_raw     <- varianceStabilizingTransformation(dds_raw_vst, blind = TRUE)
colData(vst_raw)$genotype   <- coldata$genotype
colData(vst_raw)$plant_code <- coldata$plant_code

raw_pca <- build_pca_plot(
  vst_obj = vst_raw,
  coldata = coldata,
  title   = sprintf("%s — Raw (no correction)", species_tag)
)

out_raw <- file.path(output_dir, sprintf("Raw_PCA_%s.png", species_tag))
ggsave(out_raw, plot = raw_pca, width = 20, height = 18, units = "cm",
       dpi = 320, bg = "white")
cat(sprintf("  Raw PCA saved: %s\n", out_raw))


# --- 8. edgeR model for RUVr residuals ---
cat("\n[7] Fitting edgeR model for RUVr...\n")

design <- model.matrix(~ coldata$treatment)
y      <- DGEList(counts = counts(dds_raw), group = coldata$treatment)
y      <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y      <- calcNormFactors(y, method = "upperquartile")
y      <- estimateGLMCommonDisp(y, design)
y      <- estimateGLMTagwiseDisp(y, design)
fit    <- glmFit(y, design)
res    <- residuals(fit, type = "deviance")
cat(sprintf("  Genes after edgeR filter: %d\n", nrow(y$counts)))


# --- 9. RUVr for each k ---
cat("\n[8] Running RUVr")

pca_plots <- list()

for (k in k_values) {
  cat(sprintf("\n  --- k = %d ---\n", k))
  
  ruv     <- RUVr(y$counts, rownames(y), k = k, res)
  
  dds_ruv <- DESeqDataSetFromMatrix(
    countData = ruv$normalizedCounts,
    colData   = coldata,
    design    = ~ treatment
  )
  vst_ruv <- varianceStabilizingTransformation(dds_ruv, blind = TRUE)
  colData(vst_ruv)$genotype   <- coldata$genotype
  colData(vst_ruv)$plant_code <- coldata$plant_code
  
  pca_plots[[paste0("k", k)]] <- build_pca_plot(
    vst_obj = vst_ruv,
    coldata = coldata,
    title   = sprintf("%s — RUVr k=%d", species_tag, k)
  )
  
  # Save outputs
  write.table(ruv$normalizedCounts,
              file = file.path(output_dir, sprintf("RUVr_k%d_%s_corrected_counts.tsv", k, species_tag)),
              sep = "\t", quote = FALSE, col.names = TRUE)
  
  write.table(assay(vst_ruv),
              file = file.path(output_dir, sprintf("RUVr_k%d_%s_vst.tsv", k, species_tag)),
              sep = "\t", quote = FALSE, col.names = TRUE)
  
  w_df     <- as.data.frame(ruv$W)
  w_df$run <- rownames(w_df)
  write.table(w_df,
              file = file.path(output_dir, sprintf("RUVr_k%d_%s_W_factors.tsv", k, species_tag)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat(sprintf("  RUVr k=%d outputs saved\n", k))
}


# --- 10. Save panel: Raw + k1 + k2 + k3 + k4 + k5 (3x2) ---
cat("\n[9] Saving comparison panel...\n")

k_plot_list <- lapply(k_values, function(k) pca_plots[[paste0("k", k)]])
n_plots     <- 1 + length(k_values)  # Raw PCA + one per k
n_cols      <- 3
n_rows      <- ceiling(n_plots / n_cols)

panel <- Reduce("+", c(list(raw_pca), k_plot_list)) +
  plot_annotation(
    title      = sprintf("%s — RUVr: Raw vs k=%s", species_tag, paste(k_values, collapse = ", ")),
    tag_levels = "A",
    theme = theme(
      plot.title      = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.tag        = element_text(size = 16, face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  ) +
  plot_layout(ncol = n_cols, nrow = n_rows)

panel_width  <- n_cols * 20
panel_height <- n_rows * 18

out_panel <- file.path(output_dir, sprintf("RUVr_panel_%s.png", species_tag))
ggsave(out_panel, plot = panel, width = panel_width, height = panel_height, units = "cm",
       dpi = 320, bg = "white")
cat(sprintf("  Panel saved: %s\n", out_panel))


# ============================================================
# SUMMARY TABLES
# ============================================================
cat("\n[10] Generating summary tables...\n")

# Table 1: filtering effect
all_ids_total <- rownames(y$counts)
n_genes_total <- sum(!startsWith(rownames(txi$counts), "TE_"))
n_te_total    <- sum(startsWith(rownames(txi$counts), "TE_"))
n_genes_cv    <- sum(!startsWith(rownames(counts(dds_raw)), "TE_"))
n_te_cv       <- sum(startsWith(rownames(counts(dds_raw)), "TE_"))
n_genes_kept  <- sum(!startsWith(all_ids_total, "TE_"))
n_te_kept     <- sum(startsWith(all_ids_total, "TE_"))

filter_table <- data.frame(
  Feature      = c("Protein-coding genes", "Transposable elements (TEs)", "Total"),
  Before_CV    = c(n_genes_total, n_te_total, n_genes_total + n_te_total),
  After_CV     = c(n_genes_cv,  n_te_cv,  n_genes_cv  + n_te_cv),
  After_filteredByExpression = c(n_genes_kept,  n_te_kept,  n_genes_kept  + n_te_kept),
  Retained_pct = round(100 * c(n_genes_kept, n_te_kept, n_genes_kept + n_te_kept) /
                         c(n_genes_total, n_te_total, n_genes_total + n_te_total), 1)
)
colnames(filter_table)[4] <- "Retained (%)"

out_filter <- file.path(output_dir, sprintf("summary_cv_filter_%s.tsv", species_tag))
write.table(filter_table, file = out_filter, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Filtering summary saved: %s\n", out_filter))
cat("\n  CV Filtering Effect:\n")
print(filter_table, row.names = FALSE)

# Table 2: species summary
species_summary <- coldata %>%
  summarise(
    Species        = species_tag,
    N_samples      = n(),
    N_genotypes    = n_distinct(genotype),
    N_experiments  = n_distinct(plant_code),
    N_bioprojects  = n_distinct(bioproject),
    Tissues        = paste(sort(unique(tissue)), collapse = ", "),
    Treatments     = paste(sort(unique(as.character(treatment))), collapse = ", "),
    k_tested       = paste(k_values, collapse = ", "),
    methods_run    = "RUVr"
  )

out_species <- file.path(output_dir, sprintf("summary_species_%s.tsv", species_tag))
write.table(species_summary, file = out_species, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Species summary saved: %s\n", out_species))
cat("\n  Species Summary:\n")
print(as.data.frame(species_summary), row.names = FALSE)

cat("\n================================================================\n")
cat("DONE\n")
cat("================================================================\n")
