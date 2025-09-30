# De Novo Methylation Analysis for Trio Studies
# Description: Identifies de novo methylation changes in parent-child trios using effect-size criteria

# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
# minfi (>= 1.40.0)
# IlluminaHumanMethylationEPICanno.ilm10b4.hg19 (>= 0.6.0)
# dplyr (>= 1.0.0)
# ggplot2 (>= 3.0.0)

library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(dplyr)
library(ggplot2)

# =============================================================================
# OUTPUT DIRECTORY SETUP
# =============================================================================

output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

analyze_de_novo_methylation <- function(beta_matrix_file, 
                                        metadata_file,
                                        output_dir = "results",
                                        beta_diff_threshold = 0.2,
                                        min_region_cpgs = 3) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Loading data...\n")
  
  # Load beta matrix (CpGs in rows, samples in columns)
  beta_matrix <- read.csv(beta_matrix_file, row.names = 1, check.names = FALSE)
  
  # Load metadata (must contain: Sample_ID, Family_ID, Relationship, Age)
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  # Validate required columns
  required_cols <- c("Sample_ID", "Family_ID", "Relationship", "Age")
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in metadata:", paste(missing_cols, collapse = ", ")))
  }
  
  # Identify complete trios
  trios <- identify_complete_trios(metadata)
  cat(sprintf("Found %d complete trios\n", nrow(trios)))
  
  if (nrow(trios) == 0) {
    stop("No complete trios found in the data!")
  }
  
  # Filter to only trio samples
  trio_samples <- unique(c(trios$Proband, trios$Mother, trios$Father))
  beta_matrix <- beta_matrix[, colnames(beta_matrix) %in% trio_samples, drop = FALSE]
  
  # Quality control
  cat("Performing quality control...\n")
  beta_clean <- quality_control_filtering(beta_matrix)
  
  # Remove age-associated probes
  cat("Removing age-associated probes...\n")
  age_probes <- identify_age_probes(beta_clean, metadata)
  beta_filtered <- beta_clean[!rownames(beta_clean) %in% age_probes, ]
  
  # Identify de novo methylation sites using effect-size criteria
  cat("Identifying de novo methylation sites...\n")
  de_novo_results <- find_de_novo_sites(beta_filtered, trios, beta_diff_threshold)
  
  cat(sprintf("Found %d de novo sites (Δβ > %g from both parents)\n", 
              nrow(de_novo_results), beta_diff_threshold))
  
  if (nrow(de_novo_results) == 0) {
    cat("No de novo sites found with current threshold.\n")
    return(list(results = de_novo_results, summary = data.frame()))
  }
  
  # Perform regional analysis
  cat("Performing regional analysis...\n")
  regional_results <- analyze_regions(beta_filtered, de_novo_results, trios, min_region_cpgs)
  
  # Save results
  write.csv(de_novo_results, 
            file.path(output_dir, "de_novo_methylation_sites.csv"),
            row.names = FALSE)
  
  # Generate summary statistics
  summary_stats <- generate_summary_statistics(de_novo_results, trios)
  write.csv(summary_stats,
            file.path(output_dir, "family_summary.csv"),
            row.names = FALSE)
  
  # Save regional results
  if (length(regional_results$regions) > 0) {
    for (i in seq_along(regional_results$regions)) {
      region_name <- names(regional_results$regions)[i]
      write.csv(regional_results$regions[[i]], 
                file.path(output_dir, paste0("region_", region_name, ".csv")), 
                row.names = FALSE)
    }
  }
  
  if (nrow(regional_results$gene_summary) > 0) {
    write.csv(regional_results$gene_summary, 
              file.path(output_dir, "gene_region_summary.csv"), 
              row.names = FALSE)
  }
  
  # Create visualizations
  create_summary_plots(de_novo_results, output_dir)
  
  cat("Analysis complete! Results saved to:", output_dir, "\n")
  
  return(list(
    results = de_novo_results,
    summary = summary_stats,
    regions = regional_results
  ))
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

identify_complete_trios <- function(metadata) {
  trios <- metadata %>%
    filter(Relationship %in% c("Proband", "Mother", "Father")) %>%
    group_by(Family_ID) %>%
    filter(n_distinct(Relationship) == 3) %>%
    summarise(
      Proband = Sample_ID[Relationship == "Proband"],
      Mother = Sample_ID[Relationship == "Mother"],
      Father = Sample_ID[Relationship == "Father"],
      .groups = 'drop'
    )
  return(trios)
}

quality_control_filtering <- function(beta_matrix) {
  # Remove probes with >10% missing values
  missing_threshold <- 0.1 * ncol(beta_matrix)
  keep_probes <- rowSums(is.na(beta_matrix)) <= missing_threshold
  beta_filtered <- beta_matrix[keep_probes, ]
  
  # Get annotation
  data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  # Remove SNP-affected and sex chromosome probes
  snp_probes <- rownames(anno)[!is.na(anno$Probe_SNPs) | !is.na(anno$CpG_SNPs)]
  sex_probes <- rownames(anno)[anno$chr %in% c("chrX", "chrY")]
  
  exclude_probes <- unique(c(snp_probes, sex_probes))
  beta_filtered <- beta_filtered[!rownames(beta_filtered) %in% exclude_probes, ]
  
  cat(sprintf("Retained %d/%d probes after QC\n", 
              nrow(beta_filtered), nrow(beta_matrix)))
  
  return(beta_filtered)
}

identify_age_probes <- function(beta_matrix, metadata) {
  # Load known age-associated probes from epigenetic clocks
  horvath_probes <- c("cg16867657", "cg06493994", "cg26928153", "cg03193687", "cg09809672",
                      "cg07553761", "cg12821256", "cg25809905", "cg17861230", "cg26614073",
                      "cg16419235", "cg00481951", "cg09017213", "cg12054453", "cg16054275",
                      "cg00945507", "cg17505783", "cg00864867", "cg03636183", "cg09719725",
                      "cg02085953", "cg13870906", "cg14361627", "cg16867657", "cg09809672",
                      "cg06493994", "cg07553761", "cg12821256", "cg25809905", "cg17861230",
                      "cg00481951", "cg09017213", "cg12054453", "cg16054275", "cg00945507",
                      "cg17505783", "cg00864867", "cg03636183", "cg09719725", "cg02085953",
                      "cg13870906", "cg14361627", "cg02233149", "cg25148589", "cg01968287",
                      "cg04528819", "cg09809672", "cg19761273", "cg23995914", "cg05442902",
                      "cg25148589", "cg23706838", "cg02085953", "cg00945507", "cg16867657",
                      "cg04528819", "cg12821256", "cg19761273", "cg23995914", "cg05442902",
                      "cg25148589", "cg23706838", "cg15804973", "cg22736354", "cg24768561",
                      "cg03636183", "cg00864867", "cg04875128", "cg20822990", "cg09719725",
                      "cg13870906", "cg14361627", "cg02233149", "cg01968287", "cg15804973",
                      "cg22736354", "cg24768561", "cg04875128", "cg20822990", "cg14614643",
                      "cg04528819", "cg19761273", "cg23995914", "cg05442902", "cg25148589",
                      "cg23706838", "cg15804973", "cg22736354", "cg24768561", "cg04875128",                          "cg20822990", "cg14614643", "cg06639320", "cg12054453", "cg16054275",
                      "cg17505783", "cg09017213", "cg00481951", "cg06639320", "cg16419235")
  hannum_probes <- c("cg09729957", "cg17471102", "cg09809672", "cg03193687", "cg06493994",
                     "cg16867657", "cg07553761", "cg25809905", "cg26928153", "cg12821256",
                     "cg06639320", "cg12054453", "cg16054275", "cg17505783", "cg09017213",
                     "cg00481951", "cg16419235", "cg00945507", "cg00864867", "cg03636183",
                     "cg09719725", "cg02085953", "cg13870906", "cg14361627", "cg02233149",
                     "cg25148589", "cg01968287", "cg04528819", "cg19761273", "cg23995914",
                     "cg05442902", "cg23706838", "cg15804973", "cg22736354", "cg24768561",
                     "cg04875128", "cg20822990", "cg14614643", "cg17861230", "cg26614073",
                     "cg24385618", "cg07572967", "cg00075967", "cg03193687", "cg07380416",
                     "cg02085953", "cg17471102", "cg09729957", "cg08540945", "cg26211698",
                     "cg16054275", "cg17505783", "cg06639320", "cg12054453", "cg09017213",
                     "cg00481951", "cg16419235", "cg24385618", "cg07572967", "cg00075967",
                     "cg07380416", "cg08540945", "cg26211698", "cg05442902", "cg24768561",
                     "cg15804973", "cg22736354", "cg04875128", "cg20822990", "cg14614643")
  
  age_probes <- unique(c(horvath_probes, hannum_probes))
  age_probes <- intersect(age_probes, rownames(beta_matrix))
  
  return(age_probes)
}

find_de_novo_sites <- function(beta_matrix, trios, threshold = 0.2) {
  results <- data.frame()
  
  for (i in 1:nrow(trios)) {
    family_id <- trios$Family_ID[i]
    proband <- trios$Proband[i]
    mother <- trios$Mother[i]
    father <- trios$Father[i]
    
    # Skip if samples not in beta matrix
    if (!all(c(proband, mother, father) %in% colnames(beta_matrix))) {
      next
    }
    
    # Get beta values
    proband_beta <- beta_matrix[, proband]
    mother_beta <- beta_matrix[, mother]
    father_beta <- beta_matrix[, father]
    
    # Find de novo sites (different from BOTH parents)
    diff_mother <- abs(proband_beta - mother_beta)
    diff_father <- abs(proband_beta - father_beta)
    
    de_novo_mask <- (diff_mother > threshold) & (diff_father > threshold) & 
                    !is.na(diff_mother) & !is.na(diff_father)
    
    if (sum(de_novo_mask) > 0) {
      de_novo_cpgs <- rownames(beta_matrix)[de_novo_mask]
      max_diff <- pmax(diff_mother[de_novo_mask], diff_father[de_novo_mask])
      
      family_results <- data.frame(
        Family_ID = family_id,
        CpG = de_novo_cpgs,
        Proband_Beta = proband_beta[de_novo_mask],
        Mother_Beta = mother_beta[de_novo_mask],
        Father_Beta = father_beta[de_novo_mask],
        Diff_Mother = diff_mother[de_novo_mask],
        Diff_Father = diff_father[de_novo_mask],
        Max_Diff = max_diff,
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, family_results)
    }
  }
  
  return(results)
}

analyze_regions <- function(beta_matrix, de_novo_results, trios, min_cpgs = 3) {
  if (nrow(de_novo_results) == 0) {
    return(list(regions = list(), gene_summary = data.frame()))
  }
  
  # Get genomic coordinates
  data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  anno_df <- data.frame(
    CpG_ID = rownames(anno),
    chr = as.character(anno$chr),
    pos = as.numeric(anno$pos),
    UCSC_RefGene_Name = as.character(anno$UCSC_RefGene_Name),
    UCSC_RefGene_Group = as.character(anno$UCSC_RefGene_Group),
    stringsAsFactors = FALSE
  )
  
  # Add genomic coordinates to results
  results_with_coords <- merge(de_novo_results, anno_df,
                                by.x = "CpG", by.y = "CpG_ID", all.x = TRUE)
  
  # Find regions (CpGs within 1kb)
  regions <- find_nearby_cpgs(results_with_coords, max_distance = 1000, min_cpgs = min_cpgs)
  
  # Analyze gene regions
  gene_summary <- results_with_coords %>%
    filter(!is.na(UCSC_RefGene_Name)) %>%
    group_by(UCSC_RefGene_Name, UCSC_RefGene_Group) %>%
    summarise(
      n_cpgs = n(),
      families_affected = length(unique(Family_ID)),
      mean_max_diff = mean(Max_Diff, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    filter(n_cpgs >= 3) %>%
    arrange(desc(n_cpgs))
  
  return(list(
    regions = regions,
    gene_summary = gene_summary
  ))
}

find_nearby_cpgs <- function(results, max_distance = 1000, min_cpgs = 3) {
  regions <- list()
  
  results_clean <- results[!is.na(results$chr) & !is.na(results$pos), ]
  
  if (nrow(results_clean) == 0) {
    return(regions)
  }
  
  # Group by family and chromosome
  families <- unique(results_clean$Family_ID)
  
  for (family_id in families) {
    family_data <- results_clean[results_clean$Family_ID == family_id, ]
    
    for (chr in unique(family_data$chr)) {
      chr_data <- family_data[family_data$chr == chr, ]
      chr_data <- chr_data[order(chr_data$pos), ]
      
      if (nrow(chr_data) < min_cpgs) next
      
      i <- 1
      region_id <- 1
      
      while (i <= nrow(chr_data)) {
        current_region <- chr_data[i, , drop = FALSE]
        j <- i + 1
        
        while (j <= nrow(chr_data)) {
          if (is.na(chr_data$pos[j]) || is.na(chr_data$pos[j-1])) {
            break
          }
          
          if ((chr_data$pos[j] - chr_data$pos[j-1]) <= max_distance) {
            current_region <- rbind(current_region, chr_data[j, , drop = FALSE])
            j <- j + 1
          } else {
            break
          }
        }
        
        if (nrow(current_region) >= min_cpgs) {
          region_name <- paste0("Family_", family_id, "_", chr, "_region_", region_id)
          regions[[region_name]] <- current_region
          region_id <- region_id + 1
        }
        
        i <- j
      }
    }
  }
  
  return(regions)
}

generate_summary_statistics <- function(de_novo_results, trios) {
  if (nrow(de_novo_results) == 0) {
    return(data.frame(
      Family_ID = trios$Family_ID,
      N_DeNovo_Sites = 0
    ))
  }
  
  summary_stats <- de_novo_results %>%
    group_by(Family_ID) %>%
    summarise(
      N_DeNovo_Sites = n(),
      Mean_Effect_Size = mean(Max_Diff),
      Median_Effect_Size = median(Max_Diff),
      N_Hypermethylated = sum(Proband_Beta > pmax(Mother_Beta, Father_Beta)),
      N_Hypomethylated = sum(Proband_Beta < pmin(Mother_Beta, Father_Beta)),
      .groups = 'drop'
    )
  
  # Add families with no de novo sites
  all_families <- data.frame(Family_ID = trios$Family_ID)
  summary_stats <- merge(all_families, summary_stats, by = "Family_ID", all.x = TRUE)
  summary_stats[is.na(summary_stats)] <- 0
  
  return(summary_stats)
}

create_summary_plots <- function(de_novo_results, output_dir) {
  # Sites per family
  family_counts <- de_novo_results %>%
    group_by(Family_ID) %>%
    summarise(Count = n(), .groups = 'drop')
  
  p1 <- ggplot(family_counts, aes(x = factor(Family_ID), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "De Novo Methylation Sites per Family",
         x = "Family ID", y = "Number of Sites") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "sites_per_family.png"), p1, 
         width = 8, height = 6, dpi = 300)
  
  # Effect size distribution
  p2 <- ggplot(de_novo_results, aes(x = Max_Diff)) +
    geom_histogram(bins = 30, fill = "coral", alpha = 0.7) +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    labs(title = "Distribution of Effect Sizes",
         x = "Maximum Beta Difference", y = "Count") +
    theme_minimal()
  
  ggsave(file.path(output_dir, "effect_size_distribution.png"), p2,
         width = 6, height = 4, dpi = 300)
}

# =============================================================================
# USAGE EXAMPLE
# =============================================================================

results <- analyze_de_novo_methylation(
   beta_matrix_file = "beta_matrix.csv",
   metadata_file = "sample_metadata.csv",
   output_dir = "results",
   beta_diff_threshold = 0.2,
   min_region_cpgs = 3
)
