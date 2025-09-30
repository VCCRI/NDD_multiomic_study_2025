# Generic Epigenetic Age Analysis Script
# Calculate epigenetic age using methylation data and compare groups

# =============================================================================
# REQUIRED PACKAGES
# =============================================================================
# dplyr (>= 1.0.0)
# ggplot2 (>= 3.0.0)
# readr (>= 2.0.0)
# tidyr (>= 1.0.0)
# methylclock (>= 0.7.0)
# lme4 (>= 1.1-0)
# lmerTest (>= 3.0.0)

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(methylclock)
library(lme4)
library(lmerTest)

# =============================================================================
# OUTPUT DIRECTORY SETUP
# =============================================================================

output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# DATA LOADING
# =============================================================================

# Load methylation beta values (CpGs as rows, samples as columns)
beta_matrix <- read_csv("beta_matrix.csv")
family_info <- read_csv("family_info.csv")

# Find common samples between datasets
sample_ids_beta <- colnames(beta_matrix)
individual_ids_family <- family_info$`Individual ID`
common_samples <- intersect(sample_ids_beta, individual_ids_family)

if (length(common_samples) == 0) {
  stop("No matching samples found between datasets!")
}

# Filter datasets to matching samples only
beta_matrix_filtered <- beta_matrix[, common_samples]
family_info_filtered <- family_info[family_info$`Individual ID` %in% common_samples, ]

# Check for required Family column
if (!"Family" %in% colnames(family_info_filtered)) {
  stop("Family column is required in family_info.csv!")
}

# =============================================================================
# EPIGENETIC AGE CALCULATION
# =============================================================================

# Calculate epigenetic ages using Hannum clock
hannum_age <- DNAmAge(beta_matrix_filtered,
                      clocks = "Hannum",
                      age = family_info_filtered$Age,
                      cell.count = FALSE)

# Combine results and calculate age acceleration
analysis_data <- family_info_filtered %>%
  mutate(`Individual ID` = as.character(`Individual ID`)) %>%
  left_join(data.frame(
    Individual_ID = common_samples,
    Hannum_Age = hannum_age$Hannum,
    stringsAsFactors = FALSE
  ) %>% mutate(Individual_ID = as.character(Individual_ID)),
  by = c("Individual ID" = "Individual_ID")) %>%
  mutate(
    Hannum_Acceleration = Hannum_Age - Age
  ) %>%
  filter(!is.na(Hannum_Age))

# =============================================================================
# STATISTICAL ANALYSIS
# =============================================================================

# Compare epigenetic vs chronological age
hannum_vs_chrono <- t.test(analysis_data$Hannum_Age, analysis_data$Age, paired = TRUE)

# Group comparisons (modify group variable as needed)
if ("Group" %in% colnames(analysis_data)) {
  hannum_group_test <- aov(Hannum_Acceleration ~ Group, data = analysis_data)
}

# Mixed effects models accounting for family structure
hannum_mixed <- lmer(Hannum_Acceleration ~ Group + Age + Sex + (1|Family),
                     data = analysis_data)

# =============================================================================
# VISUALIZATION
# =============================================================================

# Epigenetic age vs chronological age correlation
p1 <- ggplot(analysis_data, aes(x = Age, y = Hannum_Age)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "Hannum Epigenetic Age vs Chronological Age",
       x = "Chronological Age (years)",
       y = "Hannum Epigenetic Age (years)") +
  theme_bw()

# Age acceleration boxplots (if Group column exists)
if ("Group" %in% colnames(analysis_data)) {
  p2 <- ggplot(analysis_data, aes(x = Group, y = Hannum_Acceleration)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    labs(title = "Hannum Age Acceleration by Group",
         x = "Group",
         y = "Age Acceleration (years)") +
    theme_bw()
}

# =============================================================================
# SAVE RESULTS
# =============================================================================

write_csv(analysis_data, file.path(output_dir, "epigenetic_age_results.csv"))

ggsave(file.path(output_dir, "hannum_vs_chronological_age.png"), p1, width = 8, height = 6, dpi = 300)

if (exists("p2")) {
  ggsave(file.path(output_dir, "hannum_acceleration_by_group.png"), p2, width = 8, height = 6, dpi = 300)
}

cat("Epigenetic age analysis complete. Results saved to", output_dir, "\n")
