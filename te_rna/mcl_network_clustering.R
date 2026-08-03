library(data.table)
library(igraph)
library(ggplot2)

# ==============================================================================
# LOCAL DEFAULTS (overridden by a CLI flag)
# ==============================================================================
DEFAULTS <- list(
  edge_file            = "~/networks/edges_Sitalica_k1_r0.8_q0.05.tsv",
  output_prefix        = "~/modules/mcl_sitalica_k1",
  group_name           = "S. italica",
  inflation            = 2,    # too many tiny modules -> lower (1.4, 1.8); too few huge -> raise (2.5, 4.0)
  min_module_size      = 2,    # modules smaller than this -> "Unassigned"
  min_module_size_plot = 10,   # minimum module size shown in the plot
  cores                = 1,     # threads for C MCL (-te) and data.table
  mcl_bin              = "mcl"  # name or full path of the C MCL binary (e.g. /usr/bin/mcl)
)

# ==============================================================================
# CLI ARGUMENT PARSING
#   Example:
#     Rscript mcl_clustering.r \
#         --edge_file  /path/edges.tsv \
#         --output_prefix /path/mcl_specie_k \
#         --group_name specie \
#         --inflation  2 \
#         --cores      50
#   With no flags, it just runs the DEFAULTS block.
# ==============================================================================
parse_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit  <- which(args == flag)
  if (length(hit) && hit < length(args)) {
    return(args[hit + 1])
  }
  default   # hard-coded fallback from DEFAULTS
}

# --- I/O ----------------------------------------------------------------------
EDGE_FILE  <- parse_arg("--edge_file",     DEFAULTS$edge_file)
PREFIX     <- parse_arg("--output_prefix", DEFAULTS$output_prefix)
GROUP_NAME <- parse_arg("--group_name",    DEFAULTS$group_name)

# --- Tunable parameters -------------------------------------------------------
INFLATION            <- as.numeric(parse_arg("--inflation",            DEFAULTS$inflation))
MIN_MODULE_SIZE      <- as.integer(parse_arg("--min_module_size",      DEFAULTS$min_module_size))
MIN_MODULE_SIZE_PLOT <- as.integer(parse_arg("--min_module_size_plot", DEFAULTS$min_module_size_plot))
NUM_CORES            <- as.integer(parse_arg("--cores",                DEFAULTS$cores))
MCL_BIN              <- parse_arg("--mcl_bin", DEFAULTS$mcl_bin)

setDTthreads(NUM_CORES)

# Resolve the C MCL binary to a full path that R can definitely reach.
# Sys.which() uses R's own PATH, so this catches the case where the shell finds
# 'mcl' but R (e.g. launched from RStudio) does not.
MCL_PATH <- if (file.exists(MCL_BIN)) MCL_BIN else unname(Sys.which(MCL_BIN))
if (!nzchar(MCL_PATH)) {
  stop(sprintf(paste0(
    "C MCL binary '%s' not found on R's PATH.\n",
    "  - Confirm it is installed:            which mcl\n",
    "  - Check what R sees:                   Sys.which(\"mcl\")\n",
    "  - If the shell finds it but R doesn't, pass the full path via --mcl_bin\n",
    "    (e.g. --mcl_bin /usr/bin/mcl) or set it in the DEFAULTS block."), MCL_BIN))
}

# Basic guard: make sure the edge file actually exists before doing any work
if (!file.exists(EDGE_FILE)) {
  stop(sprintf("Edge file not found: %s\n(Check the --edge_file flag or the DEFAULTS block.)", EDGE_FILE))
}

cat("Run configuration:\n")
cat(sprintf("  edge_file            : %s\n", EDGE_FILE))
cat(sprintf("  output_prefix        : %s\n", PREFIX))
cat(sprintf("  group_name           : %s\n", GROUP_NAME))
cat(sprintf("  mcl_bin              : %s\n", MCL_PATH))
cat(sprintf("  inflation            : %.2f\n", INFLATION))
cat(sprintf("  min_module_size      : %d\n", MIN_MODULE_SIZE))
cat(sprintf("  min_module_size_plot : %d\n", MIN_MODULE_SIZE_PLOT))
cat(sprintf("  cores                : %d\n", NUM_CORES))

# --- msc_palette (custom project colours) ------------------------------------
msc_palette <- c("#ba4134", "#bf7e46", "#8894b7", "#949941",
                 "#f5c74d", "#ddc5a9", "#5d556a", "#1a657a")

# Node-type colours: distinct, colour-blind-distinguishable choices from msc_palette
# Gene = teal (#1a657a), TE = rust red (#ba4134)
NODE_TYPE_COLOURS <- c(gene = "#1a657a", TE = "#ba4134")
# Shapes as a redundant (colour-independent) encoding for scatter/point plots
NODE_TYPE_SHAPES  <- c(gene = 16, TE = 17)   # 16 = filled circle, 17 = filled triangle

# ==============================================================================
# HELPER: classify a node as TE or gene
#   Priority 1: explicit type column from the edge file (gene / TE)
#   Priority 2: fallback to the ^TE_ prefix convention
# ==============================================================================
classify_node_type <- function(node_names, type_lookup = NULL) {
  out <- rep(NA_character_, length(node_names))
  if (!is.null(type_lookup)) {
    out <- unname(type_lookup[node_names])
  }
  # Fallback / normalisation via the TE_ prefix convention
  needs_fill <- is.na(out)
  out[needs_fill] <- ifelse(grepl("^TE_", node_names[needs_fill]), "TE", "gene")
  # Normalise any stray casing/labels coming from the file
  out <- ifelse(grepl("^TE", out, ignore.case = TRUE) | grepl("^TE_", node_names),
                "TE", "gene")
  out
}

# ==============================================================================
# HELPER: Run native C MCL implementation with multi-threading
# ==============================================================================
run_mcl_multicore <- function(g, inflation, cores, mcl_path = "mcl") {
  cat(sprintf("  Executing original C MCL with %d threads (Inflation = %.2f)...\n", cores, inflation))
  
  # 1. Extract edge data frame to write out in ABC format (NodeA NodeB Weight)
  edges_df <- as_data_frame(g, what = "edges")
  
  # 2. Configure temporary files for interprocess communication
  tmp_in  <- tempfile(fileext = ".abc")
  tmp_out <- tempfile(fileext = ".mcl")
  
  # 3. Write native edge list format for MCL
  fwrite(data.table(edges_df$from, edges_df$to, edges_df$weight),
         file = tmp_in, sep = "\t", col.names = FALSE)
  
  # 4. Invoke the C-binary via system call using the multithreading flag (-te)
  status <- system2(mcl_path,
                    args = c(tmp_in, "--abc", "-I", inflation, "-te", cores, "-o", tmp_out),
                    stdout = FALSE, stderr = FALSE)
  
  if (status != 0) {
    stop("MCL execution failed. Ensure that the 'mcl' binary is installed and accessible in your system PATH.")
  }
  
  if (!file.exists(tmp_out)) {
    stop("MCL completed but output cluster file was not found.")
  }
  
  # 5. Read resulting clusters (each line contains tab-separated node names)
  lines <- readLines(tmp_out)
  
  # 6. Reconstruct an igraph compatible named membership vector
  node_names <- V(g)$name
  mem <- integer(length(node_names))
  names(mem) <- node_names
  
  for (i in seq_along(lines)) {
    if (trimws(lines[i]) == "") next
    nodes <- strsplit(lines[i], "\t")[[1]]
    mem[nodes] <- i
  }
  
  # Handle fallback edge cases: if any nodes were completely skipped by MCL,
  # assign them to single-element orphan modules
  unassigned_nodes <- names(mem)[mem == 0]
  if (length(unassigned_nodes) > 0) {
    max_cluster <- length(lines)
    for (j in seq_along(unassigned_nodes)) {
      mem[unassigned_nodes[j]] <- max_cluster + j
    }
  }
  
  # Clean up system temporary files
  unlink(c(tmp_in, tmp_out))
  
  # Calculate modularity of the MCL configuration to ensure downstream compatibility
  q <- igraph::modularity(g, mem)
  
  list(membership = mem, modularity = q)
}

# ==============================================================================
# MAIN: cluster one network
# ==============================================================================
cluster_network <- function(group_name, edge_file, prefix,
                            inflation, min_module_size, min_module_size_plot,
                            n_cores, mcl_path = "mcl") {
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("Clustering (MCL): ", group_name, "\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")
  
  # ── 1. Read edge list ─────────────────────────────────────────────
  # Real schema: Node1  Node2  weight  pvalue  qvalue  Node1_type  Node2_type
  cat("Reading edge list...\n")
  edges <- fread(edge_file, header = TRUE,
                 select = c("Node1", "Node2", "weight", "Node1_type", "Node2_type"))
  
  cat(sprintf("  Edges: %s\n",
              format(nrow(edges), big.mark = ",")))
  cat(sprintf("  Raw weight range: [%.6f, %.6f]\n",
              min(edges$weight), max(edges$weight)))
  
  # ── 1b. Preserve node-type information (gene / TE) ───────────────
  # Build a node -> type lookup from both endpoints before we drop the columns
  type_dt <- unique(rbindlist(list(
    edges[, .(node = Node1, type = Node1_type)],
    edges[, .(node = Node2, type = Node2_type)]
  )))
  # A node should have a single consistent type; keep the first occurrence
  type_dt <- type_dt[, .(type = type[1]), by = node]
  type_lookup <- setNames(type_dt$type, type_dt$node)
  
  # ── 1c. Use abs(weight) to preserve ALL strong correlations ──────
  #   Both strong positive and strong negative co-expression become strong,
  #   positive MCL edge weights (the C MCL treats weight as attraction strength).
  edges[, weight := abs(weight)]
  cat(sprintf("  |weight| range after abs(): [%.6f, %.6f]\n",
              min(edges$weight), max(edges$weight)))
  
  # ── 2. Build weighted igraph ──────────────────────────────────────
  cat("Building igraph object...\n")
  g <- graph_from_data_frame(edges[, .(Node1, Node2, weight)],
                             directed = FALSE)
  g <- simplify(g,
                remove.multiple = TRUE,
                remove.loops    = TRUE,
                edge.attr.comb  = list(weight = "max"))
  
  cat(sprintf("  Nodes: %s | Edges: %s\n",
              format(vcount(g), big.mark = ","),
              format(ecount(g), big.mark = ",")))
  rm(edges); gc()
  
  # Attach node type as a vertex attribute (gene / TE)
  V(g)$node_type <- classify_node_type(V(g)$name, type_lookup)
  n_te_total   <- sum(V(g)$node_type == "TE")
  n_gene_total <- sum(V(g)$node_type == "gene")
  cat(sprintf("  Node types -> genes: %s | TEs: %s\n",
              format(n_gene_total, big.mark = ","),
              format(n_te_total,   big.mark = ",")))
  
  # ── 3. Run MCL ────────────────────────────────────────────────────
  partition     <- run_mcl_multicore(g, inflation, n_cores, mcl_path)
  n_raw_modules <- length(unique(partition$membership))
  q             <- partition$modularity
  cat(sprintf("  Raw modules: %d | Modularity Q: %.4f\n", n_raw_modules, q))
  
  # ── 4. Build membership table ─────────────────────────────────────
  membership_dt <- data.table(
    node       = V(g)$name,
    node_type  = V(g)$node_type,
    module_raw = as.integer(partition$membership)
  )
  
  # Size-rank modules: largest = Module_001, etc.
  mod_sizes  <- membership_dt[, .N, by = module_raw][order(-N)]
  large_mods <- mod_sizes[N >= min_module_size, module_raw]
  
  rank_map <- data.table(
    module_raw  = large_mods,
    module_name = sprintf("Module_%03d", seq_along(large_mods))
  )
  
  membership_dt <- merge(membership_dt, rank_map, by = "module_raw", all.x = TRUE)
  membership_dt[is.na(module_name), module_name := "Unassigned"]
  membership_dt[, module_raw := NULL]
  
  # ── 5. Add node-level metrics ─────────────────────────────────────
  node_str <- strength(g, vids = V(g), weights = E(g)$weight)
  node_deg <- degree(g)
  
  membership_dt[, strength := node_str[node]]
  membership_dt[, degree   := node_deg[node]]
  setorder(membership_dt, module_name, -strength)
  
  # ── 6. Module summary (all modules, including small ones) ─────────
  #   Now also breaks each module down into gene vs TE counts.
  summary_dt <- membership_dt[, .(
    n_nodes = .N,
    n_genes = sum(node_type == "gene"),
    n_TEs   = sum(node_type == "TE")
  ), by = module_name][order(-n_nodes)]
  summary_dt[, modularity_Q    := round(q, 4)]
  
  # Kept column name as 'resolution' for file schema compatibility, but populated with inflation value
  summary_dt[, resolution      := inflation]
  summary_dt[, min_module_size := min_module_size]
  summary_dt[, shown_in_plot   := n_nodes >= min_module_size_plot & module_name != "Unassigned"]
  
  n_modules_all    <- sum(summary_dt$module_name != "Unassigned")
  n_modules_plot   <- sum(summary_dt$shown_in_plot)
  n_unassigned     <- membership_dt[module_name == "Unassigned", .N]
  assigned_modules <- summary_dt[module_name != "Unassigned"]
  
  cat(sprintf("  Total modules (>= %d nodes): %d\n", min_module_size, n_modules_all))
  cat(sprintf("  Modules shown in plot (>= %d nodes): %d\n", min_module_size_plot, n_modules_plot))
  cat(sprintf("  Unassigned nodes:    %d\n", n_unassigned))
  cat(sprintf("  Largest module:      %d nodes\n", max(assigned_modules$n_nodes)))
  cat(sprintf("  Median module size:  %.0f nodes\n", median(assigned_modules$n_nodes)))
  
  # ── 7. Hub nodes — top 10 per module by weighted strength ─────────
  #   node_type is carried through so hubs are labelled gene / TE.
  hubs_dt <- membership_dt[module_name != "Unassigned",
                           head(.SD, 10),
                           by      = module_name,
                           .SDcols = c("node", "node_type", "strength", "degree")]
  hubs_dt[, hub_rank := seq_len(.N), by = module_name]
  
  # ── 8. Plot data table (long format: one row per module × node_type) ──
  plot_modules <- summary_dt[shown_in_plot == TRUE][order(-n_nodes)]
  module_levels <- plot_modules$module_name
  
  plot_dt <- melt(
    plot_modules[, .(module_name, gene = n_genes, TE = n_TEs)],
    id.vars       = "module_name",
    variable.name = "node_type",
    value.name    = "n_nodes"
  )
  plot_dt[, module_name := factor(module_name, levels = module_levels)]
  plot_dt[, node_type   := factor(node_type, levels = c("gene", "TE"))]
  
  # ── 9. Save all output files ──────────────────────────────────────
  out_membership <- paste0(prefix, "_membership.tsv")
  out_summary    <- paste0(prefix, "_module_summary.tsv")
  out_hubs       <- paste0(prefix, "_hub_nodes.tsv")
  out_plot_data  <- paste0(prefix, "_plot_data.tsv")
  
  fwrite(membership_dt[, .(node, node_type, module_name, strength, degree)],
         file = out_membership, sep = "\t", quote = FALSE)
  
  fwrite(summary_dt[, .(module = module_name, n_nodes, n_genes, n_TEs,
                        modularity_Q, resolution, min_module_size, shown_in_plot)],
         file = out_summary, sep = "\t", quote = FALSE)
  
  fwrite(hubs_dt,
         file = out_hubs, sep = "\t", quote = FALSE)
  
  fwrite(plot_dt[, .(module = as.character(module_name), node_type, n_nodes)],
         file = out_plot_data, sep = "\t", quote = FALSE)
  
  cat("  Saved:", basename(out_membership), "\n")
  cat("  Saved:", basename(out_summary),    "\n")
  cat("  Saved:", basename(out_hubs),       "\n")
  cat("  Saved:", basename(out_plot_data),  "\n")
  
# ==============================================================================
# ENTRY POINT
#   One network per invocation, driven entirely by CLI arguments. To process
#   several species/networks on the cluster, submit one job per edge file
#   (e.g. an SGE array job varying --edge_file / --output_prefix / --group_name).
# ==============================================================================
result <- cluster_network(
  group_name           = GROUP_NAME,
  edge_file            = EDGE_FILE,
  prefix               = PREFIX,
  inflation            = INFLATION,
  min_module_size      = MIN_MODULE_SIZE,
  min_module_size_plot = MIN_MODULE_SIZE_PLOT,
  n_cores              = NUM_CORES,
  mcl_path             = MCL_PATH
)

rds_out <- paste0(PREFIX, "_results.rds")
saveRDS(result, file = rds_out)
cat("\nSaved RDS:", basename(rds_out), "\n")
cat("Network clustered via MCL.\n")
