################################################################################
# DMR ANALYSIS PIPELINE
################################################################################

# Required libraries
required_packages <- c(
  "minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "IlluminaHumanMethylationEPICmanifest", "limma", "DMRcate", "missMethyl",
  "maxprobes", "EpiDISH", "illuminaio", "dplyr", "ggplot2", "pheatmap",
  "RColorBrewer", "org.Hs.eg.db", "clusterProfiler", "ReactomePA", "DOSE"
)

# Check and load packages
cat("Checking required packages...\n")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall with BiocManager::install(c('", paste(missing_packages, collapse = "', '"), "'))")
}

suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(IlluminaHumanMethylationEPICmanifest)
  library(limma)
  library(DMRcate)
  library(missMethyl)
  library(maxprobes)
  library(EpiDISH)
  library(illuminaio)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(DOSE)
})

################################################################################
# CONFIGURATION
################################################################################

# Parse command line arguments or use defaults
args <- commandArgs(trailingOnly = TRUE)

# Configuration parameters
CONFIG <- list(
  # Input/output directories (can be overridden via command line)
  base_dir = ifelse(length(args) > 0, args[1], getwd()),
  sample_info_file = ifelse(length(args) > 1, args[2], "sample_information.csv"),
  output_dir = ifelse(length(args) > 2, args[3], "episignature_output"),
  
  # IDAT file locations (relative to base_dir)
  idat_dirs = list(
    GSE97362 = "geo_raw/GSE97362/idats",  # 450K array, Batch 4
    GSE74432 = "geo_raw/GSE74432/idats",  # 450K array, Batch 3
    Your_Study = "my_idats"               # EPIC array, Batch 1
  ),
  
  # Analysis parameters
  analysis_batches = c(1, 3, 4),  # Blood samples only
  qc_threshold = 0.01,            # Detection p-value threshold
  dmr_fdr_threshold = 0.05,       # FDR for DMR calling
  delta_beta_threshold = 0.2,     # Minimum beta difference for DMRs
  n_cores = 1                     # Parallel processing cores
)

# Create output directories
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "plots"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CONFIG$output_dir, "individual_profiles"), showWarnings = FALSE, recursive = TRUE)
dir.create("intermediate_data", showWarnings = FALSE, recursive = TRUE)

cat("\n=== EPISIGNATURE ANALYSIS PIPELINE ===\n")
cat("Individual patient DMR analysis (Levy et al. 2022)\n")
cat("Output directory:", CONFIG$output_dir, "\n\n")

################################################################################
# PART 1: DATA LOADING AND QC
################################################################################

################################################################################
# Load and filter sample information
################################################################################

cat("=== LOADING SAMPLE INFORMATION ===\n")

if (!file.exists(CONFIG$sample_info_file)) {
  stop("Sample information file not found: ", CONFIG$sample_info_file,
       "\n\nExpected columns: Sample, Batch, Group, Age, Sex")
}

sample_info_all <- read.csv(CONFIG$sample_info_file, stringsAsFactors = FALSE)

# Validate required columns
required_cols <- c("Sample", "Batch", "Group")
missing_cols <- setdiff(required_cols, colnames(sample_info_all))
if (length(missing_cols) > 0) {
  stop("Missing required columns in sample_information.csv: ", 
       paste(missing_cols, collapse = ", "))
}

cat("Total samples loaded:", nrow(sample_info_all), "\n")
cat("\nBatch distribution:\n")
print(table(sample_info_all$Batch))
cat("\nGroup distribution:\n")
print(table(sample_info_all$Group))

# Filter to analysis cohort
cat("\n=== FILTERING TO ANALYSIS COHORT ===\n")
sample_info_filtered <- sample_info_all[sample_info_all$Batch %in% CONFIG$analysis_batches, ]

# Define analysis groups
sample_info_filtered$Analysis_Group <- ifelse(
  sample_info_filtered$Group == "CHD_proband", "CHD_proband",
  ifelse(sample_info_filtered$Group == "NDD_proband", "NDD_proband",
         ifelse(sample_info_filtered$Group == "Control", "Control", "Exclude"))
)

# Filter to analysis groups
sample_info_filtered <- sample_info_filtered[
  sample_info_filtered$Analysis_Group != "Exclude", ]

cat("\nAnalysis cohort after filtering:\n")
print(table(sample_info_filtered$Analysis_Group, sample_info_filtered$Batch))

samples_to_keep <- sample_info_filtered$Sample

################################################################################
# Load IDAT files
################################################################################

cat("\n=== LOADING IDAT FILES ===\n")

load_filtered_idats <- function(base_dir, dataset_name, array_type, samples_to_keep) {
  idat_path <- file.path(base_dir, CONFIG$idat_dirs[[dataset_name]])
  
  if (!dir.exists(idat_path)) {
    cat("  Directory not found:", idat_path, "- skipping", dataset_name, "\n")
    return(NULL)
  }
  
  targets <- data.frame(
    Basename = list.files(idat_path, pattern = "_Grn.idat$", 
                          full.names = TRUE, recursive = TRUE)
  )
  
  if (nrow(targets) == 0) {
    cat("  No IDAT files found in", idat_path, "- skipping", dataset_name, "\n")
    return(NULL)
  }
  
  targets$Basename <- gsub("_Grn.idat$", "", targets$Basename)
  targets$Sample_Name <- basename(targets$Basename)
  targets$Array <- array_type
  targets$Dataset <- dataset_name
  
  targets <- targets[targets$Sample_Name %in% samples_to_keep, ]
  
  if (nrow(targets) == 0) {
    cat("  No matching samples for", dataset_name, "\n")
    return(NULL)
  }
  
  cat("  Loading", nrow(targets), "samples from", dataset_name, "(", array_type, ")\n")
  RGset <- read.metharray.exp(targets = targets, force = TRUE)
  return(list(RGset = RGset, targets = targets))
}

# Load datasets
result_97362 <- load_filtered_idats(CONFIG$base_dir, "GSE97362", "450K", samples_to_keep)
result_74432 <- load_filtered_idats(CONFIG$base_dir, "GSE74432", "450K", samples_to_keep)
result_epic <- load_filtered_idats(CONFIG$base_dir, "Your_Study", "EPIC", samples_to_keep)

################################################################################
# Quality control
################################################################################

cat("\n=== QUALITY CONTROL ===\n")

qc_and_filter_samples <- function(result, dataset_name, threshold = CONFIG$qc_threshold) {
  if (is.null(result)) return(NULL)
  
  RGset <- result$RGset
  targets <- result$targets
  
  # Calculate detection p-values
  detP <- detectionP(RGset)
  mean_detP <- colMeans(detP)
  failed_samples <- mean_detP > threshold
  
  cat("\n", dataset_name, ":\n", sep = "")
  cat("  Mean detP range:", 
      sprintf("%.4f - %.4f", min(mean_detP), max(mean_detP)), "\n")
  cat("  Failed samples (detP >", threshold, "):", sum(failed_samples), "\n")
  
  # Extract negative control probe intensities
  cat("  Extracting control probe intensities...\n")
  ctrl_probes <- getProbeInfo(RGset, type = "Control")
  negative_ctrl_idx <- which(ctrl_probes$Type == "NEGATIVE")
  
  ctrl_intensities_log <- NULL
  if (length(negative_ctrl_idx) > 0) {
    red_ctrl <- getRed(RGset)[ctrl_probes$Address[negative_ctrl_idx], , drop = FALSE]
    green_ctrl <- getGreen(RGset)[ctrl_probes$Address[negative_ctrl_idx], , drop = FALSE]
    ctrl_intensities <- (red_ctrl + green_ctrl) / 2
    ctrl_intensities_log <- log2(ctrl_intensities + 1)
    cat("    Found", length(negative_ctrl_idx), "negative control probes\n")
  }
  
  # Filter samples
  keep <- !failed_samples
  RGset_filtered <- RGset[, keep]
  detP_filtered <- detP[, keep]
  targets_filtered <- targets[keep, ]
  
  if (!is.null(ctrl_intensities_log)) {
    ctrl_intensities_log <- ctrl_intensities_log[, keep, drop = FALSE]
  }
  
  cat("  Samples retained:", ncol(RGset_filtered), "\n")
  
  return(list(
    RGset = RGset_filtered,
    detP = detP_filtered,
    targets = targets_filtered,
    control_probes = ctrl_intensities_log
  ))
}

result_97362 <- qc_and_filter_samples(result_97362, "GSE97362")
result_74432 <- qc_and_filter_samples(result_74432, "GSE74432")
result_epic <- qc_and_filter_samples(result_epic, "Your_Study")

################################################################################
# Normalization (preprocessNoob for cross-platform compatibility)
################################################################################

cat("\n=== NORMALIZATION ===\n")
cat("Using preprocessNoob for optimal cross-platform analysis\n\n")

normalize_dataset <- function(result, dataset_name) {
  if (is.null(result)) return(NULL)
  
  cat("Normalizing", dataset_name, "...\n")
  mSet <- preprocessNoob(result$RGset, dyeCorr = TRUE, verbose = FALSE)
  cat("  Complete:", ncol(mSet), "samples\n")
  
  return(list(
    mSet = mSet,
    detP = result$detP,
    targets = result$targets,
    control_probes = result$control_probes
  ))
}

result_97362 <- normalize_dataset(result_97362, "GSE97362")
result_74432 <- normalize_dataset(result_74432, "GSE74432")
result_epic <- normalize_dataset(result_epic, "Your_Study")

################################################################################
# Merge datasets and filter probes
################################################################################

cat("\n=== MERGING DATASETS ===\n")

# Combine datasets (keeping only common probes for cross-platform compatibility)
datasets <- list(result_97362, result_74432, result_epic)
datasets <- datasets[!sapply(datasets, is.null)]

if (length(datasets) == 0) {
  stop("No datasets loaded successfully!")
}

# Find common probes
common_probes <- Reduce(intersect, lapply(datasets, function(x) rownames(x$mSet)))
cat("Common probes across platforms:", length(common_probes), "\n")

# Extract beta values and merge
all_betas <- do.call(cbind, lapply(datasets, function(x) {
  getBeta(x$mSet)[common_probes, ]
}))

all_targets <- do.call(rbind, lapply(datasets, function(x) x$targets))
all_detP <- do.call(cbind, lapply(datasets, function(x) {
  x$detP[common_probes, ]
}))

cat("Total samples after merging:", ncol(all_betas), "\n")

# Filter probes based on detection p-value
failed_probes <- rowMeans(all_detP > CONFIG$qc_threshold) > 0.1
cat("Probes failing detection (>10% samples with detP >", CONFIG$qc_threshold, "):", 
    sum(failed_probes), "\n")

all_betas_clean <- all_betas[!failed_probes, ]
cat("Probes retained after filtering:", nrow(all_betas_clean), "\n")

# Remove problematic probes (cross-reactive, SNPs, sex chromosomes)
cat("\n=== FILTERING PROBLEMATIC PROBES ===\n")

# Get annotation
if (any(grepl("cg", rownames(all_betas_clean)))) {
  anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
} else {
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
}

# Remove sex chromosome probes
sex_probes <- rownames(anno)[anno$chr %in% c("chrX", "chrY")]
all_betas_clean <- all_betas_clean[!rownames(all_betas_clean) %in% sex_probes, ]
cat("After removing sex chromosome probes:", nrow(all_betas_clean), "\n")

# Remove SNP probes
snp_probes <- rownames(anno)[!is.na(anno$Probe_rs) | !is.na(anno$CpG_rs)]
all_betas_clean <- all_betas_clean[!rownames(all_betas_clean) %in% snp_probes, ]
cat("After removing SNP-affected probes:", nrow(all_betas_clean), "\n")

# Remove cross-reactive probes (maxprobes package)
if (requireNamespace("maxprobes", quietly = TRUE)) {
  if (any(grepl("cg", rownames(all_betas_clean)))) {
    xreactive <- maxprobes::xreactive_probes(array_type = "450K")
  } else {
    xreactive <- maxprobes::xreactive_probes(array_type = "EPIC")
  }
  all_betas_clean <- all_betas_clean[!rownames(all_betas_clean) %in% xreactive, ]
  cat("After removing cross-reactive probes:", nrow(all_betas_clean), "\n")
}

################################################################################
# Batch correction using ComBat
################################################################################

cat("\n=== BATCH CORRECTION ===\n")

# Merge sample info with targets
all_targets_merged <- merge(all_targets, sample_info_filtered, 
                            by.x = "Sample_Name", by.y = "Sample", 
                            all.x = TRUE)

# Ensure same order
all_targets_merged <- all_targets_merged[match(colnames(all_betas_clean), 
                                               all_targets_merged$Sample_Name), ]

# Apply ComBat
if ("Batch" %in% colnames(all_targets_merged)) {
  batch_vector <- as.factor(all_targets_merged$Batch)
  
  # Create model matrix with biological variables
  mod <- model.matrix(~ Analysis_Group, data = all_targets_merged)
  
  cat("Applying ComBat batch correction...\n")
  cat("  Batches:", paste(unique(batch_vector), collapse = ", "), "\n")
  
  all_betas_corrected <- limma::removeBatchEffect(
    all_betas_clean,
    batch = batch_vector,
    design = mod
  )
  
  # Ensure beta values remain in [0, 1] range
  all_betas_corrected[all_betas_corrected < 0] <- 0
  all_betas_corrected[all_betas_corrected > 1] <- 1
  
  cat("Batch correction complete\n")
} else {
  cat("No batch information found, skipping batch correction\n")
  all_betas_corrected <- all_betas_clean
}

# Save intermediate data
save(all_betas_corrected, all_targets_merged, 
     file = "intermediate_data/processed_data.RData")
cat("\nIntermediate data saved: intermediate_data/processed_data.RData\n")

################################################################################
# PART 2: INDIVIDUAL PATIENT DMR ANALYSIS
################################################################################

cat("\n=== INDIVIDUAL PATIENT DMR ANALYSIS ===\n")
cat("Approach: Levy et al. 2022 - individual vs matched controls\n\n")

# Get control samples
control_samples <- all_targets_merged$Sample_Name[
  all_targets_merged$Analysis_Group == "Control"
]
cat("Control samples available:", length(control_samples), "\n")

# Identify patient samples
patient_samples_CHD <- all_targets_merged$Sample_Name[
  all_targets_merged$Analysis_Group == "CHD_proband"
]
patient_samples_NDD <- all_targets_merged$Sample_Name[
  all_targets_merged$Analysis_Group == "NDD_proband"
]

cat("CHD patient samples:", length(patient_samples_CHD), "\n")
cat("NDD patient samples:", length(patient_samples_NDD), "\n\n")

################################################################################
# Function: Individual patient DMR analysis
################################################################################

analyze_individual_dmrs <- function(patient_id, patient_betas, control_betas,
                                   fdr_threshold = CONFIG$dmr_fdr_threshold,
                                   delta_beta_threshold = CONFIG$delta_beta_threshold) {
  
  cat("Analyzing patient:", patient_id, "\n")
  
  # Combine patient with controls
  combined_betas <- cbind(patient_betas, control_betas)
  design <- model.matrix(~ factor(c(1, rep(0, ncol(control_betas)))))
  colnames(design) <- c("Intercept", "Patient")
  
  # Limma differential methylation
  fit <- lmFit(combined_betas, design)
  fit <- eBayes(fit)
  
  # Get results
  results <- topTable(fit, coef = "Patient", number = Inf, sort.by = "none")
  results$probe <- rownames(results)
  results$delta_beta <- rowMeans(patient_betas) - rowMeans(control_betas)
  
  # Filter by FDR and delta beta
  sig_probes <- results[results$adj.P.Val < fdr_threshold & 
                       abs(results$delta_beta) > delta_beta_threshold, ]
  
  cat("  Significant CpG sites (FDR <", fdr_threshold, ", |Δβ| >", 
      delta_beta_threshold, "):", nrow(sig_probes), "\n")
  
  if (nrow(sig_probes) == 0) {
    cat("  No significant probes found\n")
    return(NULL)
  }
  
  # Annotate probes
  anno_sub <- anno[match(sig_probes$probe, rownames(anno)), ]
  sig_probes$chr <- anno_sub$chr
  sig_probes$pos <- anno_sub$pos
  sig_probes$gene <- anno_sub$UCSC_RefGene_Name
  sig_probes$feature <- anno_sub$UCSC_RefGene_Group
  
  # Call DMRs using DMRcate
  cat("  Calling DMRs with DMRcate...\n")
  
  # Prepare annotation
  myAnnotation <- cpg.annotate(
    datatype = "array",
    object = combined_betas,
    what = "Beta",
    arraytype = ifelse(any(grepl("cg", rownames(combined_betas))), "450K", "EPIC"),
    analysis.type = "differential",
    design = design,
    coef = "Patient",
    fdr = fdr_threshold
  )
  
  # Find DMRs
  dmrs <- tryCatch({
    DMRcate::dmrcate(myAnnotation, lambda = 1000, C = 2, min.cpgs = 3)
  }, error = function(e) {
    cat("    DMRcate failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(dmrs) || is.null(dmrs$results) || nrow(dmrs$results) == 0) {
    cat("  No DMRs identified\n")
    return(list(probes = sig_probes, dmrs = NULL))
  }
  
  dmr_results <- extractRanges(dmrs, genome = "hg19")
  cat("  DMRs identified:", length(dmr_results), "\n")
  
  # Convert to data frame
  dmr_df <- as.data.frame(dmr_results)
  dmr_df$patient_id <- patient_id
  
  return(list(
    probes = sig_probes,
    dmrs = dmr_df,
    n_probes = nrow(sig_probes),
    n_dmrs = nrow(dmr_df)
  ))
}

################################################################################
# Run analysis for all patients
################################################################################

# Analyze CHD patients
cat("\n=== ANALYZING CHD PATIENTS ===\n")
CHD_results <- list()

for (patient_id in patient_samples_CHD) {
  patient_betas <- all_betas_corrected[, patient_id, drop = FALSE]
  control_betas <- all_betas_corrected[, control_samples]
  
  result <- analyze_individual_dmrs(patient_id, patient_betas, control_betas)
  
  if (!is.null(result)) {
    CHD_results[[patient_id]] <- result
    
    # Save individual results
    write.csv(result$probes, 
              file.path(CONFIG$output_dir, "individual_profiles",
                       paste0(patient_id, "_significant_probes.csv")),
              row.names = FALSE)
    
    if (!is.null(result$dmrs)) {
      write.csv(result$dmrs,
                file.path(CONFIG$output_dir, "individual_profiles",
                         paste0(patient_id, "_DMRs.csv")),
                row.names = FALSE)
    }
  }
}

# Analyze NDD patients
cat("\n=== ANALYZING NDD PATIENTS ===\n")
NDD_results <- list()

for (patient_id in patient_samples_NDD) {
  patient_betas <- all_betas_corrected[, patient_id, drop = FALSE]
  control_betas <- all_betas_corrected[, control_samples]
  
  result <- analyze_individual_dmrs(patient_id, patient_betas, control_betas)
  
  if (!is.null(result)) {
    NDD_results[[patient_id]] <- result
    
    # Save individual results
    write.csv(result$probes,
              file.path(CONFIG$output_dir, "individual_profiles",
                       paste0(patient_id, "_significant_probes.csv")),
              row.names = FALSE)
    
    if (!is.null(result$dmrs)) {
      write.csv(result$dmrs,
                file.path(CONFIG$output_dir, "individual_profiles",
                         paste0(patient_id, "_DMRs.csv")),
                row.names = FALSE)
    }
  }
}

################################################################################
# Aggregate results across patients
################################################################################

cat("\n=== AGGREGATING RESULTS ===\n")

# Combine all DMRs
all_CHD_dmrs <- do.call(rbind, lapply(CHD_results, function(x) x$dmrs))
all_NDD_dmrs <- do.call(rbind, lapply(NDD_results, function(x) x$dmrs))

if (!is.null(all_CHD_dmrs)) {
  write.csv(all_CHD_dmrs,
            file.path(CONFIG$output_dir, "CHD_all_patient_DMRs.csv"),
            row.names = FALSE)
  cat("CHD: Total DMRs across all patients:", nrow(all_CHD_dmrs), "\n")
}

if (!is.null(all_NDD_dmrs)) {
  write.csv(all_NDD_dmrs,
            file.path(CONFIG$output_dir, "NDD_all_patient_DMRs.csv"),
            row.names = FALSE)
  cat("NDD: Total DMRs across all patients:", nrow(all_NDD_dmrs), "\n")
}

