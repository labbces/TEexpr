# ==============================================================================
# MCL post-processing & TE-focused visualisation
# ==============================================================================
library(data.table)
library(ggplot2)

# ==============================================================================
# LOCAL DEFAULTS  (overridden by CLI flags)
# ==============================================================================
DEFAULTS <- list(
  membership_file = "/home/cunha/results/networks/modules/mcl_sitalica_k1_membership.tsv",
  edge_file       = "/home/cunha/results/networks/edges_Sitalica_k1_r0.8_q0.05.tsv",  # ORIGINAL signed edges
  output_prefix   = "/home/cunha/results/networks/modules/mcl_sitalica_k1",
  group_name      = "S. italica",
  enrich_alpha    = 0.05,  # BH q-value threshold to flag a module as TE-enriched (table only)
  min_hub_module_size = 10,  # ignore modules smaller than this when picking the "most central TE"
  cores           = 1
)

# ==============================================================================
# CLI ARGUMENT PARSING (parse_arg convention with DEFAULTS fallback)
#   Example:
#     Rscript mcl_postprocess.r \
#         --membership_file /path/mcl_sitalica_k1_membership.tsv \
#         --edge_file       /path/edges_Sitalica_k1_r0.8_q0.05.tsv \
#         --output_prefix   /path/mcl_sitalica_k1 \
#         --group_name      "S. italica" \
#         --cores 4
# ==============================================================================
parse_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit  <- which(args == flag)
  if (length(hit) && hit < length(args)) return(args[hit + 1])
  default
}

MEMBERSHIP_FILE <- parse_arg("--membership_file", DEFAULTS$membership_file)
EDGE_FILE       <- parse_arg("--edge_file",       DEFAULTS$edge_file)
PREFIX          <- parse_arg("--output_prefix",   DEFAULTS$output_prefix)
GROUP_NAME      <- parse_arg("--group_name",      DEFAULTS$group_name)
ENRICH_ALPHA    <- as.numeric(parse_arg("--enrich_alpha", DEFAULTS$enrich_alpha))
MIN_HUB_MODULE_SIZE <- as.integer(parse_arg("--min_hub_module_size", DEFAULTS$min_hub_module_size))
NUM_CORES       <- as.integer(parse_arg("--cores", DEFAULTS$cores))

setDTthreads(NUM_CORES)

for (f in c(MEMBERSHIP_FILE, EDGE_FILE)) {
  if (!file.exists(f)) stop(sprintf("Input file not found: %s", f))
}

cat("Post-processing configuration:\n")
cat(sprintf("  membership_file : %s\n", MEMBERSHIP_FILE))
cat(sprintf("  edge_file       : %s\n", EDGE_FILE))
cat(sprintf("  output_prefix   : %s\n", PREFIX))
cat(sprintf("  group_name      : %s\n", GROUP_NAME))
cat(sprintf("  enrich_alpha    : %.3f\n", ENRICH_ALPHA))
cat(sprintf("  min_hub_module_size : %d\n", MIN_HUB_MODULE_SIZE))

# --- msc_palette (custom project colours) ------------------------------------
msc_palette <- c("#ba4134", "#bf7e46", "#8894b7", "#949941",
                 "#f5c74d", "#ddc5a9", "#5d556a", "#1a657a")

# Node-type colours + redundant shapes (colour-independent encoding)
NODE_TYPE_COLOURS <- c(gene = "#1a657a", TE = "#ba4134")   # gene=teal, TE=rust
NODE_TYPE_SHAPES  <- c(gene = 16, TE = 17)                 # circle / triangle

# Shared theme so every figure matches dissertation formatting ----------------
theme_msc <- function(base = 13) {
  theme_classic(base_size = base, base_family = "Helvetica") +
    theme(
      plot.title    = element_text(face = "bold", size = 26, hjust = 0.5, family = "Helvetica",
                                   margin = margin(b = 10)),
      plot.subtitle = element_text(size = 13, hjust = 0.5, family = "Helvetica", margin = margin(b = 8)),
      plot.margin   = margin(t = 15, r = 20, b = 10, l = 10),
      axis.title    = element_text(size = 16, family = "Helvetica"),
      axis.text     = element_text(size = 12, family = "Helvetica"),
      legend.title  = element_text(size = 16, family = "Helvetica"),
      legend.text   = element_text(size = 14, family = "Helvetica"),
      legend.position = "top",
      panel.grid.major = element_line(colour = "grey87", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

# Helper: normalise a type label / node name to "gene" or "TE"
norm_type <- function(type_vec, node_vec) {
  ifelse(grepl("^TE", type_vec, ignore.case = TRUE) | grepl("^TE_", node_vec), "TE", "gene")
}

# ==============================================================================
# 1. READ INPUTS
# ==============================================================================
cat("\nReading membership table...\n")
mem <- fread(MEMBERSHIP_FILE, header = TRUE)   # node, node_type, module_name, strength, degree
mem[, node_type := norm_type(node_type, node)]

cat("Reading signed edge list...\n")
edges <- fread(EDGE_FILE, header = TRUE,
               select = c("Node1", "Node2", "weight", "Node1_type", "Node2_type"))

# Global counts (population for the enrichment test = ALL nodes in the network)
N_total <- nrow(mem)
M_TE    <- mem[node_type == "TE",   .N]
N_gene  <- mem[node_type == "gene", .N]
p_TE    <- M_TE / N_total
cat(sprintf("  Nodes: %s | genes: %s | TEs: %s | global TE fraction: %.4f\n",
            format(N_total, big.mark = ","), format(N_gene, big.mark = ","),
            format(M_TE, big.mark = ","), p_TE))

# ==============================================================================
# 2. PER-MODULE MASTER TABLE
# ==============================================================================
cat("\nBuilding per-module master table...\n")

# 2a. Composition + connectivity summaries per module.
#     as.numeric() keeps column types consistent across groups (empty median -> NA_real_).
mod_tab <- mem[, .(
    n_nodes              = .N,
    n_genes              = sum(node_type == "gene"),
    n_TEs                = sum(node_type == "TE"),
    median_strength_gene = round(as.numeric(median(strength[node_type == "gene"])), 3),
    median_strength_TE   = round(as.numeric(median(strength[node_type == "TE"])),   3),
    median_degree_gene   = round(as.numeric(median(degree[node_type == "gene"])),   1),
    median_degree_TE     = round(as.numeric(median(degree[node_type == "TE"])),     1)
  ), by = module_name]

mod_tab[, prop_TE := round(n_TEs / n_nodes, 4)]

# 2b. TE enrichment via hypergeometric test (over-representation), real modules only.
real_mods <- mod_tab[module_name != "Unassigned"]
real_mods[, expected_TE     := round(n_nodes * p_TE, 2)]
real_mods[, fold_enrichment := round((n_TEs / n_nodes) / p_TE, 3)]
real_mods[, enrich_pvalue   := phyper(n_TEs - 1, m = M_TE, n = N_gene,
                                       k = n_nodes, lower.tail = FALSE)]
real_mods[, enrich_qvalue   := p.adjust(enrich_pvalue, method = "BH")]
real_mods[, TE_enriched     := enrich_qvalue < ENRICH_ALPHA & fold_enrichment > 1]

enrich_cols <- c("module_name", "expected_TE", "fold_enrichment",
                 "enrich_pvalue", "enrich_qvalue", "TE_enriched")
mod_tab <- merge(mod_tab, real_mods[, ..enrich_cols], by = "module_name", all.x = TRUE)

# 2c. Drop the "Unassigned" bucket — it is NOT a real module (nodes MCL left
#     ungrouped / below MIN_MODULE_SIZE). They remain in the *_membership.tsv.
n_unassigned_nodes <- mod_tab[module_name == "Unassigned", n_nodes]
mod_tab <- mod_tab[module_name != "Unassigned"]
if (length(n_unassigned_nodes)) {
  cat(sprintf("  Dropped 'Unassigned' bucket from table (%d ungrouped nodes; kept in membership file).\n",
              n_unassigned_nodes))
}

setorder(mod_tab, -n_nodes)
out_table <- paste0(PREFIX, "_module_TE_table.tsv")
fwrite(mod_tab, out_table, sep = "\t", quote = FALSE)
cat("  Saved:", basename(out_table), "\n")

n_enriched <- sum(mod_tab$TE_enriched, na.rm = TRUE)
cat(sprintf("  TE-enriched modules (q < %.3f, fold > 1): %d\n", ENRICH_ALPHA, n_enriched))

# ==============================================================================
# 3. FIGURE 1 — module ranking: number of TEs and TE proportion
# ==============================================================================
cat("\nFigure 1: module ranking by size, TE count, and TE proportion...\n")

N_TOP <- 15
rank_dt <- mod_tab[n_TEs >= 1][order(-n_TEs)][seq_len(min(N_TOP, .N))]
module_order <- rev(rank_dt$module_name)   # reversed so largest ends up on top after coord_flip
rank_dt[, module_name := factor(module_name, levels = module_order)]

# Built manually (not via melt) so each panel can carry its own display label:
# module size is log10-transformed for the bar height (huge range across
# modules) but labelled with the real node count so it stays readable.
panel_size  <- rank_dt[, .(module_name, panel = "Module size (log10 nodes)",
                           value = log10(n_nodes), label = format(n_nodes, big.mark = ","))]
panel_count <- rank_dt[, .(module_name, panel = "TE count",
                           value = n_TEs, label = as.character(n_TEs))]
panel_prop  <- rank_dt[, .(module_name, panel = "TE proportion",
                           value = prop_TE, label = sprintf("%.1f%%", prop_TE * 100))]
rank_long <- rbindlist(list(panel_size, panel_count, panel_prop))
rank_long[, panel := factor(panel, levels = c("Module size (log10 nodes)", "TE count", "TE proportion"))]

PANEL_COLOURS <- c("Module size (log10 nodes)" = "#8894b7", "TE count" = "#ba4134", "TE proportion" = "#bf7e46")

p_rank <- ggplot(rank_long, aes(x = module_name, y = value, fill = panel)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.15, alpha = 0.9) +
  geom_text(aes(label = label), hjust = -0.15, size = 3, family = "Helvetica") +
  coord_flip(clip = "off") +
  facet_wrap(~ panel, scales = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.22))) +
  scale_fill_manual(values = PANEL_COLOURS, guide = "none") +
  labs(
    title    = paste("Top", nrow(rank_dt), "modules by TE count \u2014", GROUP_NAME),
    x = "Module", y = NULL
  ) +
  theme_msc() +
  theme(strip.text = element_text(size = 13, face = "bold", family = "Helvetica"),
        strip.background = element_rect(fill = "#ddc5a9", colour = NA),
        axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(1.3, "lines"))

ggsave(paste0(PREFIX, "_module_TE_ranking.pdf"), p_rank, width = 12, height = 6.5)
ggsave(paste0(PREFIX, "_module_TE_ranking.png"), p_rank, width = 12, height = 6.5, dpi = 300)
cat("  Saved:", basename(paste0(PREFIX, "_module_TE_ranking.pdf")), "\n")

# ==============================================================================
# 4. FIGURE 2 — module composition (stacked bars): gene + TE counts per module
# ==============================================================================
cat("Figure 2: module composition (stacked bars, TE-containing modules only)...\n")

comp_dt <- mod_tab[n_TEs >= 1][order(-n_nodes)]
comp_dt[, module_name := factor(module_name, levels = module_name)]
n_te_containing <- nrow(comp_dt)
cat(sprintf("  TE-containing modules shown: %d (of %d total)\n", n_te_containing, nrow(mod_tab)))

comp_long <- melt(comp_dt[, .(module_name, Gene = n_genes, TE = n_TEs)],
                  id.vars = "module_name", variable.name = "node_type", value.name = "n")
comp_long[, node_type := factor(node_type, levels = c("Gene", "TE"))]

sizes_sorted <- sort(comp_dt$n_nodes, decreasing = TRUE)
largest <- sizes_sorted[1]
second  <- if (length(sizes_sorted) >= 2) sizes_sorted[2] else largest

axis_mode <- "linear"

# NOTE: comp_long's node_type levels are "Gene"/"TE" (from the melt() column
# names below), so the fill palette must use matching (capitalised) names
STACK_COLOURS <- c(Gene = "#1a657a", TE = "#ba4134")

p_comp <- ggplot(comp_long, aes(x = module_name, y = n, fill = node_type)) +
  geom_col(alpha = 0.9, colour = "grey20", linewidth = 0.15) +
  scale_fill_manual(values = STACK_COLOURS, breaks = c("Gene", "TE")) +
  labs(
    title = sprintf("Module composition\n(TE-containing modules)"),
    subtitle = GROUP_NAME,
    x = "Module (sorted by size)", y = "Number of nodes"
  ) +
  theme_msc() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

use_break <- largest > 2 * second
if (use_break && requireNamespace("ggbreak", quietly = TRUE)) {
  # Wider break: start just above the second-largest bar and jump almost all
  # the way to the top, so the compressed region spans roughly "~400 to ~3000"
  # instead of a more conservative break.
  brk_low  <- ceiling(second * 1.05)
  brk_high <- floor(largest * 0.97)
  if (brk_high > brk_low) {
    p_comp <- p_comp + ggbreak::scale_y_break(c(brk_low, brk_high), scales = 0.4, space = 0.4, symbol = "slash")
    axis_mode <- "broken"
    cat(sprintf("  Axis break applied between %d and %d (linear scale, proportions preserved).\n",
                brk_low, brk_high))
  }
}

axis_lbl <- switch(axis_mode,
                   broken = "broken y-axis (linear)",
                   log10  = "log10 y-axis (install 'ggbreak' for a true break)",
                   "linear y-axis")
p_comp <- p_comp + labs(subtitle = "Modules containing transposable elements")

ggsave(paste0(PREFIX, "_module_composition.pdf"), p_comp, width = 10, height = 7.5)
ggsave(paste0(PREFIX, "_module_composition.png"), p_comp, width = 10, height = 7.5, dpi = 300)
cat("  Saved:", basename(paste0(PREFIX, "_module_composition.pdf")), "\n")

# ==============================================================================
# 5. FIGURE 3 — most central TE per module (one TE per TE-containing module)
#    For every module with >= 1 TE AND at least --min_hub_module_size nodes,
#    pick the single TE node with the highest weighted strength (its most
#    central TE). Small modules are excluded here: with few nodes, a TE can
#    look "central" just because there is little competition, not because it
#    is genuinely well-connected - centrality there is not yet meaningful.
# ==============================================================================
cat("Figure 3: most central TE per module...\n")
cat(sprintf("  Restricting to modules with >= %d nodes (--min_hub_module_size)\n", MIN_HUB_MODULE_SIZE))

eligible_modules <- mod_tab[n_nodes >= MIN_HUB_MODULE_SIZE, module_name]
n_excluded <- mod_tab[n_TEs >= 1 & n_nodes < MIN_HUB_MODULE_SIZE, .N]
if (n_excluded > 0) {
  cat(sprintf("  Excluded %d small TE-containing module(s) below the size threshold.\n", n_excluded))
}

te_nodes <- mem[node_type == "TE" & module_name %in% eligible_modules]
top_te_per_module <- te_nodes[order(-strength), .SD[1], by = module_name]
top_te_per_module <- top_te_per_module[order(-strength)]

n_te_modules_total <- nrow(top_te_per_module)
TOP_N_HUB <- 20
top_te_per_module <- top_te_per_module[seq_len(min(TOP_N_HUB, .N))]
n_te_modules <- nrow(top_te_per_module)
cat(sprintf("  Modules containing >= 1 TE (eligible): %d | showing top %d by weighted strength\n",
            n_te_modules_total, n_te_modules))

if (n_te_modules == 0) {
  cat("  No TE-containing modules found - skipping Figure 3.\n")
} else {
  hub_dt <- copy(top_te_per_module)
  hub_dt[, module_lbl := factor(module_name, levels = rev(module_name))]  # largest strength on top

  p_hub <- ggplot(hub_dt, aes(x = strength, y = module_lbl)) +
    geom_segment(aes(x = 0, xend = strength, y = module_lbl, yend = module_lbl),
                colour = "#8894b7", linewidth = 0.6) +
    geom_point(size = 3.2, colour = "#ba4134", shape = 17) +
    geom_text(aes(label = node), hjust = -0.15, size = 3, family = "Helvetica") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
    labs(
      title    = paste("Most central TE per module \u2014", GROUP_NAME),
      subtitle = sprintf("Top %d modules by weighted strength (of %d eligible, >= %d nodes)",
                         n_te_modules, n_te_modules_total, MIN_HUB_MODULE_SIZE),
      x = "Weighted strength (sum of |r|)", y = "Module"
    ) +
    theme_msc() +
    theme(axis.text.y = element_text(size = 9))

  fig_height <- max(4, 0.28 * n_te_modules + 1.5)
  ggsave(paste0(PREFIX, "_TE_hub_per_module.pdf"), p_hub, width = 9, height = fig_height)
  ggsave(paste0(PREFIX, "_TE_hub_per_module.png"), p_hub, width = 9, height = fig_height, dpi = 300)
  cat("  Saved:", basename(paste0(PREFIX, "_TE_hub_per_module.pdf")), "\n")

  out_hub_table <- paste0(PREFIX, "_TE_hub_per_module.tsv")
  fwrite(top_te_per_module[, .(module_name, TE_node = node, strength, degree)],
         out_hub_table, sep = "\t", quote = FALSE)
  cat("  Saved:", basename(out_hub_table), "\n")
}

cat("\nAll post-processing outputs written with prefix:\n  ", PREFIX, "\n", sep = "")
cat("Done.\n")
