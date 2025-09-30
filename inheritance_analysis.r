# Inheritance Analysis
# Analyze methylation inheritance patterns in trio families with quality filtering

# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
# ggplot2 (>= 3.0.0)
# dplyr (>= 1.0.0)

library(ggplot2)
library(dplyr)

# =============================================================================
# OUTPUT DIRECTORY SETUP
# =============================================================================

output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# INHERITANCE CALCULATION FUNCTIONS
# =============================================================================

calculate_inheritance <- function(beta_matrix, family_info, target_family_ids,
                                         min_parent_diff = 0.2, require_between_parents = TRUE) {

  cat("=== INHERITANCE ANALYSIS ===\n")
  cat("Filtering criteria:\n")
  cat("1. Parents must differ by >", min_parent_diff, "\n")
  cat("2. Proband must fall between parents:", require_between_parents, "\n\n")

  results <- list()

  for (fam_id in target_family_ids) {
    cat("Processing family:", fam_id, "\n")

    fam_data <- family_info[family_info$Family_ID == fam_id, ]

    if (nrow(fam_data) != 3) {
      cat("  Warning: Family does not have exactly 3 members\n")
      next
    }

    # Get sample IDs
    proband_id <- fam_data$Sample_ID[fam_data$Relationship == "Proband"]
    mother_id <- fam_data$Sample_ID[fam_data$Relationship == "Mother"]
    father_id <- fam_data$Sample_ID[fam_data$Relationship == "Father"]

    if (!all(c(proband_id, mother_id, father_id) %in% colnames(beta_matrix))) {
      cat("  Warning: Not all family members found in beta matrix\n")
      next
    }

    # Get methylation values for all CpGs
    proband_beta <- beta_matrix[, proband_id]
    mother_beta <- beta_matrix[, mother_id]
    father_beta <- beta_matrix[, father_id]

    # Apply filtering criteria
    filtered_results <- apply_inheritance_filters(
      proband_beta, mother_beta, father_beta,
      min_parent_diff, require_between_parents
    )

    if (filtered_results$n_valid == 0) {
      cat("  Warning: No CpGs passed filtering criteria\n")
      next
    }

    # Store results
    results[[as.character(fam_id)]] <- list(
      Family_ID = fam_id,
      cpg_ids = rownames(beta_matrix)[filtered_results$valid_mask],
      n_total_cpgs = length(proband_beta),
      n_valid_cpgs = filtered_results$n_valid,
      prop_valid = filtered_results$n_valid / length(proband_beta),
      prop_closer_mother = mean(filtered_results$closer_to_mother),
      prop_closer_father = mean(filtered_results$closer_to_father),
      mean_parent_diff = mean(filtered_results$parent_diff_filtered),
      mean_inheritance_score = mean(filtered_results$inheritance_scores)
    )

    cat("  Total CpGs:", length(proband_beta), "\n")
    cat("  Valid CpGs after filtering:", filtered_results$n_valid,
        "(", round(filtered_results$n_valid/length(proband_beta)*100, 1), "%)\n")
    cat("  Maternal inheritance:", round(mean(filtered_results$closer_to_mother)*100, 1), "%\n")
    cat("  Paternal inheritance:", round(mean(filtered_results$closer_to_father)*100, 1), "%\n\n")
  }

  return(results)
}

apply_inheritance_filters <- function(proband_beta, mother_beta, father_beta,
                                      min_parent_diff, require_between_parents) {

  # Remove missing values
  complete_cases <- !is.na(proband_beta) & !is.na(mother_beta) & !is.na(father_beta)

  proband_clean <- proband_beta[complete_cases]
  mother_clean <- mother_beta[complete_cases]
  father_clean <- father_beta[complete_cases]

  # Calculate parent differences
  parent_diff <- abs(mother_clean - father_clean)

  # Filter 1: Parents must differ by minimum threshold
  sufficient_diff <- parent_diff > min_parent_diff

  # Filter 2: Proband must fall between parents (if required)
  if (require_between_parents) {
    min_parent <- pmin(mother_clean, father_clean)
    max_parent <- pmax(mother_clean, father_clean)
    between_parents <- (proband_clean >= min_parent) & (proband_clean <= max_parent)
  } else {
    between_parents <- rep(TRUE, length(proband_clean))
  }

  # Combined filter
  valid_mask_clean <- sufficient_diff & between_parents

  # Create full-length mask
  valid_mask <- rep(FALSE, length(proband_beta))
  valid_mask[complete_cases] <- valid_mask_clean

  if (sum(valid_mask_clean) == 0) {
    return(list(n_valid = 0, valid_mask = valid_mask))
  }

  # Apply filters to get final dataset
  proband_filtered <- proband_clean[valid_mask_clean]
  mother_filtered <- mother_clean[valid_mask_clean]
  father_filtered <- father_clean[valid_mask_clean]
  parent_diff_filtered <- parent_diff[valid_mask_clean]

  # Calculate distances and inheritance for filtered data
  dist_to_mother <- abs(proband_filtered - mother_filtered)
  dist_to_father <- abs(proband_filtered - father_filtered)

  closer_to_mother <- as.numeric(dist_to_mother < dist_to_father)
  closer_to_father <- as.numeric(dist_to_father < dist_to_mother)

  # Calculate inheritance scores (positive = more like mother, negative = more like father)
  inheritance_scores <- (dist_to_father - dist_to_mother) / parent_diff_filtered

  return(list(
    n_valid = sum(valid_mask_clean),
    valid_mask = valid_mask,
    parent_diff_filtered = parent_diff_filtered,
    closer_to_mother = closer_to_mother,
    closer_to_father = closer_to_father,
    inheritance_scores = inheritance_scores
  ))
}

perform_inheritance_tests <- function(results) {
  cat("=== STATISTICAL TESTS ===\n")

  test_results <- list()

  for (fam_id in names(results)) {
    fam_data <- results[[fam_id]]

    if (fam_data$n_valid_cpgs < 10) {
      cat("Family", fam_id, ": Insufficient data for testing\n")
      next
    }

    # Binomial test for inheritance bias
    maternal_count <- round(fam_data$prop_closer_mother * fam_data$n_valid_cpgs)
    test <- binom.test(x = maternal_count, n = fam_data$n_valid_cpgs, p = 0.5)

    test_results[[fam_id]] <- list(
      Family_ID = fam_data$Family_ID,
      maternal_count = maternal_count,
      total_cpgs = fam_data$n_valid_cpgs,
      prop_maternal = fam_data$prop_closer_mother,
      p_value = test$p.value,
      significant = test$p.value < 0.05,
      bias_direction = ifelse(fam_data$prop_closer_mother > 0.5, "Maternal", "Paternal")
    )

    cat("Family", fam_id, ":",
        ifelse(fam_data$prop_closer_mother > 0.5, "Maternal", "Paternal"),
        "bias, p =", round(test$p.value, 4),
        ifelse(test$p.value < 0.05, "(Significant)", "(Not significant)"), "\n")
  }

  return(test_results)
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# =============================================================================

run_inheritance_analysis <- function(beta_matrix_file, 
                                          metadata_file,
                                          target_family_ids = NULL,
                                          min_parent_diff = 0.2,
                                          output_dir = "results") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("=== INHERITANCE ANALYSIS ===\n")

  # Load data
  cat("Loading data...\n")
  beta_matrix <- read.csv(beta_matrix_file, row.names = 1, check.names = FALSE)
  family_info <- read.csv(metadata_file, stringsAsFactors = FALSE)

  # If no target families specified, use all families
  if (is.null(target_family_ids)) {
    target_family_ids <- unique(family_info$Family_ID)
  }

  cat("Target families:", paste(target_family_ids, collapse = ", "), "\n\n")

  # Run inheritance analysis
  cat("Running inheritance calculation...\n")
  results <- calculate_inheritance(
    beta_matrix, family_info, target_family_ids,
    min_parent_diff, require_between_parents = TRUE
  )

  # Perform statistical tests
  test_results <- perform_inheritance_tests(results)

  # Create summary data frame
  summary_df <- data.frame(
    Family_ID = sapply(results, function(x) x$Family_ID),
    total_cpgs = sapply(results, function(x) x$n_total_cpgs),
    valid_cpgs = sapply(results, function(x) x$n_valid_cpgs),
    prop_valid = sapply(results, function(x) x$prop_valid),
    prop_maternal = sapply(results, function(x) x$prop_closer_mother),
    mean_inheritance_score = sapply(results, function(x) x$mean_inheritance_score),
    stringsAsFactors = FALSE
  )

  # Add test results if available
  if (length(test_results) > 0) {
    test_df <- data.frame(
      Family_ID = sapply(test_results, function(x) x$Family_ID),
      p_value = sapply(test_results, function(x) x$p_value),
      significant = sapply(test_results, function(x) x$significant),
      bias_direction = sapply(test_results, function(x) x$bias_direction),
      stringsAsFactors = FALSE
    )
    summary_df <- merge(summary_df, test_df, by = "Family_ID", all.x = TRUE)
  }

  # Save results
  write.csv(summary_df, file.path(output_dir, "inheritance_summary.csv"), row.names = FALSE)

  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Results saved to:", output_dir, "\n")

  return(list(
    results = results,
    test_results = test_results,
    summary = summary_df
  ))
}

# =============================================================================
# USAGE EXAMPLE
# =============================================================================

results <- run_inheritance_analysis(
  beta_matrix_file = "beta_matrix.csv",
  metadata_file = "sample_metadata.csv",
  target_family_ids = c(1, 2, 3),  # or NULL for all families
  min_parent_diff = 0.2,
  output_dir = "results"
)