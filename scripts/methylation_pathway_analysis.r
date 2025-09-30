# Generic DNA Methylation Pathway Analysis
# Identify DMRs and perform GO enrichment analysis

library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(clusterProfiler)
library(org.Hs.eg.db)

# =============================================================================
# DATA LOADING
# =============================================================================

# Load differentially methylated region (DMR) data
# Expected columns: seqnames, start, end, Gene.Name, probeID, Annotation, direction, Entrez.ID, Gene.Type
dmr_data <- read.delim("dmr_results.txt", header = TRUE)

# =============================================================================
# COLLAPSE PROBES INTO DMRS
# =============================================================================

collapse_to_dmrs <- function(dmr_data, distance_threshold = 500) {
  # Create a GRanges object
  gr <- GRanges(
    seqnames = dmr_data$seqnames,
    ranges = IRanges(start = dmr_data$start, end = dmr_data$end),
    gene = dmr_data$Gene.Name,
    probe_id = dmr_data$probeID,
    annotation = dmr_data$Annotation,
    direction = dmr_data$direction,
    entrezID = dmr_data$Entrez.ID,
    gene_type = dmr_data$Gene.Type
  )

  # Sort and reduce GRanges by proximity
  gr <- sort(gr)
  reduced_gr <- reduce(gr, min.gapwidth = distance_threshold, with.revmap = TRUE)

  # Extract consolidated information
  dmrs <- data.frame(
    chromosome = as.character(seqnames(reduced_gr)),
    start = start(reduced_gr),
    end = end(reduced_gr),
    gene = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(unique(mcols(gr)$gene[idx]), collapse = ";")
    }),
    probes = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(mcols(gr)$probe_id[idx], collapse = ",")
    }),
    annotation = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(unique(mcols(gr)$annotation[idx]), collapse = ";")
    }),
    direction = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(unique(mcols(gr)$direction[idx]), collapse = ";")
    }),
    entrezID = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(unique(mcols(gr)$entrezID[idx]), collapse = ";")
    }),
    gene_type = sapply(mcols(reduced_gr)$revmap, function(idx) {
      paste(unique(mcols(gr)$gene_type[idx]), collapse = ";")
    })
  )

  return(dmrs)
}

# Collapse CpG probes into DMRs
dmr_results <- collapse_to_dmrs(dmr_data, distance_threshold = 500)

# Filter to keep only DMRs with multiple probes (≥5 recommended)
dmr_results_filtered <- dmr_results %>%
  mutate(probes_count = sapply(strsplit(probes, ","), length)) %>%
  filter(probes_count >= 5)

cat("Total DMRs identified:", nrow(dmr_results), "\n")
cat("DMRs with ≥5 probes:", nrow(dmr_results_filtered), "\n")

# =============================================================================
# GENE ONTOLOGY ENRICHMENT ANALYSIS
# =============================================================================

# Extract unique Entrez IDs for enrichment analysis
entrez_ids <- unique(unlist(strsplit(dmr_results_filtered$entrezID, ";")))
entrez_ids <- entrez_ids[!is.na(entrez_ids) & entrez_ids != ""]

if (length(entrez_ids) > 0) {
  # Biological Process GO enrichment
  go_bp_results <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )

  # Molecular Function GO enrichment
  go_mf_results <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )

  # Cellular Component GO enrichment
  go_cc_results <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05
  )

  cat("GO BP terms identified:", nrow(go_bp_results@result), "\n")
  cat("GO MF terms identified:", nrow(go_mf_results@result), "\n")
  cat("GO CC terms identified:", nrow(go_cc_results@result), "\n")

} else {
  cat("Warning: No valid Entrez IDs found for enrichment analysis\n")
}

# =============================================================================
# VISUALIZATION
# =============================================================================

if (exists("go_bp_results") && nrow(go_bp_results@result) > 0) {
  # Create dot plot for top GO terms
  p1 <- dotplot(go_bp_results, showCategory = 20) +
    ggtitle("GO Biological Process Enrichment")

  # Create bar plot
  p2 <- barplot(go_bp_results, showCategory = 15) +
    ggtitle("Top GO BP Terms")

  # Save plots
  ggsave("go_bp_dotplot.png", p1, width = 10, height = 8, dpi = 300)
  ggsave("go_bp_barplot.png", p2, width = 10, height = 6, dpi = 300)
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

# Save DMR results
write.csv(dmr_results_filtered, "dmr_results_filtered.csv", row.names = FALSE)

# Save GO results
if (exists("go_bp_results")) {
  write.csv(go_bp_results@result, "go_bp_enrichment.csv", row.names = FALSE)
}
if (exists("go_mf_results")) {
  write.csv(go_mf_results@result, "go_mf_enrichment.csv", row.names = FALSE)
}
if (exists("go_cc_results")) {
  write.csv(go_cc_results@result, "go_cc_enrichment.csv", row.names = FALSE)
}

cat("Analysis complete. Results saved to CSV files.\n")
