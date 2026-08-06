library(topGO)

# ============================================================================
# DEFAULTS - Edit here for local testing only.
# ============================================================================
DEFAULTS <- list(
  gene2go   = "~/terms_GO/Sbi_go_simpl_genes.txt",
  universe  = "~/terms_GO/Sbi_gene_ids.txt",
  modules   = "~/terms_GO/mcl_sbicolor_membership.tsv",
  outdir    = "~/terms_GO/topGO_results",
  prefix    = "sbi",
  ontology  = "BP",
  nodesize  = 4,
  padj      = 0.05
)

args <- commandArgs(trailingOnly = TRUE)

# reads "--flag value" from the command line, falling back to the DEFAULTS entry
parse_arg <- function(flag, default) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) stop("Missing value for ", flag)
  value <- args[i + 1]
  # keep the type of the default (numeric flags stay numeric)
  if (is.numeric(default)) as.numeric(value) else value
}

GENE2GO_FILE <- parse_arg("--gene2go",  DEFAULTS$gene2go)
UNIVERSE_FILE<- parse_arg("--universe", DEFAULTS$universe)
MODULES_FILE <- parse_arg("--modules",  DEFAULTS$modules)
OUTDIR       <- parse_arg("--outdir",   DEFAULTS$outdir)
PREFIX       <- parse_arg("--prefix",   DEFAULTS$prefix)
ONTOLOGY     <- parse_arg("--ontology", DEFAULTS$ontology)
NODESIZE     <- parse_arg("--nodesize", DEFAULTS$nodesize)
PADJ_CUTOFF  <- parse_arg("--padj",     DEFAULTS$padj)

cat("gene2GO file :", GENE2GO_FILE,  "\n")
cat("Universe file:", UNIVERSE_FILE, "\n")
cat("Modules file :", MODULES_FILE,  "\n")
cat("Output dir   :", OUTDIR,        "\n")
cat("Ontology     :", ONTOLOGY,      "\n")
cat("Node size    :", NODESIZE,      "\n")
cat("padj cutoff  :", PADJ_CUTOFF,   "\n")

# gene-to-GO mapping: two columns, gene ID and GO term
GTOGO<-read.delim2(GENE2GO_FILE, header = FALSE)

# prep table format for topGO: named list, one character vector of GO terms per gene
geneID2GO<- by(GTOGO$V2,
               GTOGO$V1,
               function(x) as.character(x))

# define universe: genes from the whole network
allGenes <- read.delim2(UNIVERSE_FILE, header = TRUE)
allGenes$GeneID <- as.factor(allGenes$GeneID)

# MCL membership table: one row per node, with its module assignment
modules <- read.delim2(MODULES_FILE, header = TRUE)

# every module holding at least one gene node (TE-only modules are skipped)
moduleNames <- sort(unique(modules$module_name[modules$node_type == "gene"]))

cat("Modules to analyse:", length(moduleNames), "\n")

# collects the number of enriched terms per module, written out at the end
summaryTable <- data.frame()

for (mod in moduleNames) {
  
  cat("\n== Processing", mod, "==\n")
  
  # genes belonging to the current module (the foreground set)
  module_ids <- as.data.frame(modules[modules$node_type == "gene" & modules$module_name == mod, "node"])
  colnames(module_ids) <- c("GeneID")
  
  # binary vector over the universe: 1 = gene in this module, 0 = background
  geneList <- factor(as.integer(allGenes$GeneID %in% module_ids$GeneID))
  names(geneList)<-allGenes$GeneID
  
  # a module with no gene matching the universe yields a single-level factor,
  # which topGO cannot test
  if (length(levels(geneList)) < 2) {
    cat("  No module gene found in the universe. Module skipped.\n")
    next
  }
  
  # topGO data preparation for the selected ontology
  module_topGO<-new("topGOdata",
                    description = mod, ontology = ONTOLOGY,
                    allGenes = geneList,
                    nodeSize = NODESIZE,
                    annot =  annFUN.gene2GO, gene2GO = geneID2GO)
  
  # classic: each GO term tested independently
  # weight01: accounts for the GO graph structure, down-weighting terms whose
  # signal comes from an already significant child term
  resultFisher  <- runTest(module_topGO, algorithm = "classic",  statistic = "fisher")
  resultWeight01<- runTest(module_topGO, algorithm = "weight01", statistic = "fisher")
  
  # term-level counts plus the raw classic p-value
  resultTable<-cbind(termStat(module_topGO),
                     score(resultFisher))
  
  colnames(resultTable)<-c('Annotated','Significant','Expected','pvalue')
  
  # multiple testing correction over the classic p-values
  resultTable$padj<-p.adjust(resultTable$pvalue, method = 'fdr')
  
  # match weight01 scores by GO ID (rownames), since their order may differ
  # from termStat; positional binding would silently misalign the columns
  w01 <- score(resultWeight01)
  resultTable$weight01 <- w01[rownames(resultTable)]
  
  # term identifiers, human-readable description and module label
  resultTable$GO.ID <- rownames(resultTable)
  resultTable$Term  <- Term(resultTable$GO.ID)
  resultTable$Module<- mod
  
  # sort by adjusted p-value and fix the column order for the output file
  resultTable <- resultTable[order(resultTable$padj, resultTable$pvalue),
                             c('Module','GO.ID','Term','Annotated','Significant',
                               'Expected','pvalue','padj','weight01')]
  
  # count enriched GO terms
  nSig <- nrow(resultTable[which(resultTable$padj < PADJ_CUTOFF),])
  cat("  Enriched GO terms (padj <", PADJ_CUTOFF, "):", nSig, "\n")
  
  # one TSV per module, holding every tested term (filter downstream if needed)
  outFile <- file.path(OUTDIR, paste0(PREFIX, "_", mod, "_enriched_go_terms.tsv"))
  write.table(resultTable, file = outFile,
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  summaryTable <- rbind(summaryTable,
                        data.frame(Module = mod,
                                   Genes_in_module = nrow(module_ids),
                                   GO_terms_tested = nrow(resultTable),
                                   Enriched = nSig))
}

# overview across modules, useful to pick which ones deserve a figure
write.table(summaryTable, file = file.path(OUTDIR, paste0(PREFIX, "_topGO_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Results written to", OUTDIR, "\n")
