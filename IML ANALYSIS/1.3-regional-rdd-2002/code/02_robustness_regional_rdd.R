# ============================================================================
# 02_robustness_regional_rdd.R
# ============================================================================
# Purpose: Robustness checks for regional RDD analysis
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# Tests:
#   1. Stability across different cutoff search ranges
#   2. Leave-one-region-out analysis
#   3. Regional heterogeneity tests
#   4. Weighted vs unweighted pooled estimates
#   5. Sensitivity to small regions
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rdrobust)
library(gridExtra)

# ============================================================================
# SETUP
# ============================================================================

# Load data
merged_entry_iml_rg <- read.csv("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

cat("========================================\n")
cat("REGIONAL RDD ROBUSTNESS CHECKS\n")
cat("========================================\n\n")

# Function to find regional cutoffs (baseline)
find_regional_cutoffs <- function(data, post_year = 2002, 
                                  cutoff_min = -1.5, cutoff_max = -0.5, 
                                  step = 0.05) {
  regions <- sort(unique(data$region41))
  regional_cutoffs <- data.frame()
  
  for (region_num in regions) {
    region_data <- data %>% 
      filter(region41 == region_num) %>%
      filter(!is.na(iml_95) & !is.na(entry_year_30))
    
    if (nrow(region_data) < 30) next
    
    region_data$post_year <- ifelse(region_data$entry_year_30 >= post_year, 1, 0)
    
    cutoff_range <- seq(cutoff_min, cutoff_max, by = step)
    max_effect <- 0
    optimal_cutoff <- NA
    
    for (cutoff in cutoff_range) {
      region_data$treatment <- ifelse(region_data$iml_95 >= cutoff, 1, 0)
      tryCatch({
        model <- lm(post_year ~ treatment, data = region_data)
        coef_val <- abs(coef(model)["treatment"])
        if (!is.na(coef_val) && coef_val > max_effect) {
          max_effect <- coef_val
          optimal_cutoff <- cutoff
        }
      }, error = function(e) NULL)
    }
    
    if (!is.na(optimal_cutoff)) {
      regional_cutoffs <- rbind(regional_cutoffs,
                                data.frame(region = region_num,
                                           cutoff = optimal_cutoff,
                                           effect = max_effect))
    }
  }
  return(regional_cutoffs)
}

# Get baseline cutoffs
baseline_cutoffs <- find_regional_cutoffs(merged_entry_iml_rg)
cat("Baseline: Found cutoffs for", nrow(baseline_cutoffs), "regions\n\n")

# ============================================================================
# TEST 1: STABILITY OF CUTOFF SELECTION
# ============================================================================

cat("------------------------------------\n")
cat("TEST 1: CUTOFF SEARCH SENSITIVITY\n")
cat("------------------------------------\n")
cat("Testing different search ranges\n\n")

# Test different search ranges
search_ranges <- list(
  "Narrow" = c(-1.3, -0.7),
  "Baseline" = c(-1.5, -0.5),
  "Wide" = c(-1.7, -0.3),
  "Shifted_down" = c(-1.6, -0.6),
  "Shifted_up" = c(-1.4, -0.4)
)

cutoff_sensitivity <- data.frame()

for (range_name in names(search_ranges)) {
  range_vals <- search_ranges[[range_name]]
  cutoffs <- find_regional_cutoffs(merged_entry_iml_rg, 
                                   cutoff_min = range_vals[1],
                                   cutoff_max = range_vals[2])
  
  # Normalize and run national RDD
  normalized_data <- merged_entry_iml_rg %>%
    mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
    inner_join(cutoffs, by = c("region41" = "region")) %>%
    mutate(normalized_iml = iml_95 - cutoff) %>%
    filter(!is.na(normalized_iml),
           entry_year_30 >= 1995 & entry_year_30 <= 2010)
  
  rd_result <- rdrobust(
    y = normalized_data$post_2002,
    x = normalized_data$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  cutoff_sensitivity <- rbind(cutoff_sensitivity, data.frame(
    range = range_name,
    min = range_vals[1],
    max = range_vals[2],
    n_regions = nrow(cutoffs),
    estimate = rd_result$coef[1],
    se = rd_result$se[1],
    pvalue = rd_result$pv[1]
  ))
  
  cat(sprintf("%-15s [%.1f, %.1f]: Est = %.3f (SE = %.3f), %d regions\n",
              range_name, range_vals[1], range_vals[2], 
              rd_result$coef[1], rd_result$se[1], nrow(cutoffs)))
}

# ============================================================================
# TEST 2: LEAVE-ONE-REGION-OUT
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 2: LEAVE-ONE-REGION-OUT\n")
cat("------------------------------------\n")
cat("Testing influence of individual regions\n\n")

# Run analysis leaving out each region
leave_out_results <- data.frame()

for (exclude_region in unique(baseline_cutoffs$region)) {
  # Get cutoffs excluding this region
  cutoffs_subset <- baseline_cutoffs %>% filter(region != exclude_region)
  
  # Normalize and run RDD
  normalized_data <- merged_entry_iml_rg %>%
    filter(region41 != exclude_region) %>%
    mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
    inner_join(cutoffs_subset, by = c("region41" = "region")) %>%
    mutate(normalized_iml = iml_95 - cutoff) %>%
    filter(!is.na(normalized_iml),
           entry_year_30 >= 1995 & entry_year_30 <= 2010)
  
  rd_result <- tryCatch({
    rdrobust(
      y = normalized_data$post_2002,
      x = normalized_data$normalized_iml,
      c = 0,
      kernel = "triangular",
      bwselect = "mserd"
    )
  }, error = function(e) NULL)
  
  if (!is.null(rd_result)) {
    leave_out_results <- rbind(leave_out_results, data.frame(
      excluded_region = exclude_region,
      estimate = rd_result$coef[1],
      se = rd_result$se[1]
    ))
  }
}

# Calculate statistics
if (nrow(leave_out_results) > 0) {
  cat("Leave-one-out estimates:\n")
  cat("  Mean:", round(mean(leave_out_results$estimate), 3), "\n")
  cat("  SD:", round(sd(leave_out_results$estimate), 3), "\n")
  cat("  Min:", round(min(leave_out_results$estimate), 3), "\n")
  cat("  Max:", round(max(leave_out_results$estimate), 3), "\n")
  
  # Plot distribution
  p_loo <- ggplot(leave_out_results, aes(x = estimate)) +
    geom_histogram(bins = 20, fill = "lightblue", color = "darkblue") +
    geom_vline(xintercept = mean(leave_out_results$estimate), 
               linetype = "dashed", color = "red") +
    labs(title = "Leave-One-Region-Out: Distribution of Estimates",
         x = "RD Estimate",
         y = "Frequency",
         caption = "Red line = mean estimate") +
    theme_minimal()
  
  ggsave("output/robustness_leave_one_out_regional.pdf", p_loo, width = 10, height = 6)
}

# ============================================================================
# TEST 3: REGIONAL HETEROGENEITY
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 3: REGIONAL HETEROGENEITY\n")
cat("------------------------------------\n\n")

# Group regions by characteristics and test differences
# First, let's group by sample size quartiles
regional_data_sizes <- merged_entry_iml_rg %>%
  group_by(region41) %>%
  summarise(n_obs = n(), .groups = 'drop')

quartiles <- quantile(regional_data_sizes$n_obs, probs = c(0.25, 0.5, 0.75))
regional_data_sizes <- regional_data_sizes %>%
  mutate(size_group = case_when(
    n_obs <= quartiles[1] ~ "Q1 (Small)",
    n_obs <= quartiles[2] ~ "Q2",
    n_obs <= quartiles[3] ~ "Q3",
    TRUE ~ "Q4 (Large)"
  ))

# Run RDD by size group
size_group_results <- data.frame()

for (group in unique(regional_data_sizes$size_group)) {
  regions_in_group <- regional_data_sizes %>%
    filter(size_group == group) %>%
    pull(region41)
  
  # Filter cutoffs and data
  cutoffs_group <- baseline_cutoffs %>%
    filter(region %in% regions_in_group)
  
  normalized_data <- merged_entry_iml_rg %>%
    filter(region41 %in% regions_in_group) %>%
    mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
    inner_join(cutoffs_group, by = c("region41" = "region")) %>%
    mutate(normalized_iml = iml_95 - cutoff) %>%
    filter(!is.na(normalized_iml),
           entry_year_30 >= 1995 & entry_year_30 <= 2010)
  
  rd_result <- rdrobust(
    y = normalized_data$post_2002,
    x = normalized_data$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  size_group_results <- rbind(size_group_results, data.frame(
    group = group,
    n_regions = length(regions_in_group),
    estimate = rd_result$coef[1],
    se = rd_result$se[1],
    pvalue = rd_result$pv[1]
  ))
  
  cat(sprintf("%-12s (%d regions): Est = %.3f (SE = %.3f)\n",
              group, length(regions_in_group), 
              rd_result$coef[1], rd_result$se[1]))
}

# ============================================================================
# TEST 4: WEIGHTED VS UNWEIGHTED
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 4: WEIGHTED POOLED ESTIMATES\n")
cat("------------------------------------\n\n")

# Create weights based on region size
weight_data <- merged_entry_iml_rg %>%
  group_by(region41) %>%
  mutate(region_n = n()) %>%
  ungroup() %>%
  mutate(weight = region_n / n())

# Normalize with baseline cutoffs
weighted_data <- weight_data %>%
  mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
  inner_join(baseline_cutoffs, by = c("region41" = "region")) %>%
  mutate(normalized_iml = iml_95 - cutoff) %>%
  filter(!is.na(normalized_iml),
         entry_year_30 >= 1995 & entry_year_30 <= 2010)

# Run weighted RDD
rd_weighted <- rdrobust(
  y = weighted_data$post_2002,
  x = weighted_data$normalized_iml,
  c = 0,
  kernel = "triangular",
  bwselect = "mserd",
  weights = weighted_data$weight
)

# Run unweighted for comparison
rd_unweighted <- rdrobust(
  y = weighted_data$post_2002,
  x = weighted_data$normalized_iml,
  c = 0,
  kernel = "triangular",
  bwselect = "mserd"
)

cat("Unweighted estimate:", round(rd_unweighted$coef[1], 3), 
    "(SE =", round(rd_unweighted$se[1], 3), ")\n")
cat("Weighted estimate:  ", round(rd_weighted$coef[1], 3), 
    "(SE =", round(rd_weighted$se[1], 3), ")\n")

# ============================================================================
# TEST 5: EXCLUDING SMALL REGIONS
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 5: SENSITIVITY TO SMALL REGIONS\n")
cat("------------------------------------\n\n")

min_obs_thresholds <- c(30, 50, 100, 200, 500)
small_region_results <- data.frame()

for (min_obs in min_obs_thresholds) {
  # Filter regions by size
  large_regions <- regional_data_sizes %>%
    filter(n_obs >= min_obs) %>%
    pull(region41)
  
  cutoffs_filtered <- baseline_cutoffs %>%
    filter(region %in% large_regions)
  
  normalized_data <- merged_entry_iml_rg %>%
    filter(region41 %in% large_regions) %>%
    mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
    inner_join(cutoffs_filtered, by = c("region41" = "region")) %>%
    mutate(normalized_iml = iml_95 - cutoff) %>%
    filter(!is.na(normalized_iml),
           entry_year_30 >= 1995 & entry_year_30 <= 2010)
  
  if (nrow(normalized_data) > 100) {
    rd_result <- rdrobust(
      y = normalized_data$post_2002,
      x = normalized_data$normalized_iml,
      c = 0,
      kernel = "triangular",
      bwselect = "mserd"
    )
    
    small_region_results <- rbind(small_region_results, data.frame(
      min_obs = min_obs,
      n_regions = length(large_regions),
      estimate = rd_result$coef[1],
      se = rd_result$se[1],
      pvalue = rd_result$pv[1]
    ))
    
    cat(sprintf("Min obs ≥ %3d: %2d regions, Est = %.3f (SE = %.3f)\n",
                min_obs, length(large_regions), 
                rd_result$coef[1], rd_result$se[1]))
  }
}

# ============================================================================
# SUMMARY PLOTS
# ============================================================================

cat("\n------------------------------------\n")
cat("CREATING SUMMARY PLOTS\n")
cat("------------------------------------\n\n")

# Plot cutoff sensitivity
if (nrow(cutoff_sensitivity) > 0) {
  p_cutoff <- ggplot(cutoff_sensitivity, aes(x = range, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), 
                  width = 0.2) +
    labs(title = "Sensitivity to Cutoff Search Range",
         x = "Search Range",
         y = "RD Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("output/robustness_cutoff_sensitivity_regional.pdf", 
         p_cutoff, width = 10, height = 6)
}

# Plot small region sensitivity
if (nrow(small_region_results) > 0) {
  p_small <- ggplot(small_region_results, aes(x = min_obs, y = estimate)) +
    geom_line(color = "blue") +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), 
                alpha = 0.2, fill = "blue") +
    labs(title = "Sensitivity to Excluding Small Regions",
         x = "Minimum Observations per Region",
         y = "RD Estimate") +
    theme_minimal()
  
  ggsave("output/robustness_small_regions_regional.pdf", 
         p_small, width = 10, height = 6)
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n========================================\n")
cat("REGIONAL ROBUSTNESS SUMMARY\n")
cat("========================================\n\n")

# Compile key results
summary_regional <- data.frame(
  Test = c("Baseline",
           "Wide search range",
           "Leave-one-out (mean)",
           "Weighted by region size",
           "Exclude regions < 100 obs"),
  Estimate = c(cutoff_sensitivity$estimate[cutoff_sensitivity$range == "Baseline"],
               cutoff_sensitivity$estimate[cutoff_sensitivity$range == "Wide"],
               mean(leave_out_results$estimate),
               rd_weighted$coef[1],
               small_region_results$estimate[small_region_results$min_obs == 100]),
  SE = c(cutoff_sensitivity$se[cutoff_sensitivity$range == "Baseline"],
         cutoff_sensitivity$se[cutoff_sensitivity$range == "Wide"],
         mean(leave_out_results$se),
         rd_weighted$se[1],
         small_region_results$se[small_region_results$min_obs == 100])
)

summary_regional$Estimate <- round(summary_regional$Estimate, 3)
summary_regional$SE <- round(summary_regional$SE, 3)

print(summary_regional)

# Save results
save(cutoff_sensitivity, leave_out_results, size_group_results, 
     small_region_results, summary_regional,
     file = "output/regional_robustness_results.RData")

write.csv(summary_regional, "output/regional_robustness_summary.csv", row.names = FALSE)

cat("\n✓ All regional robustness checks complete\n")
cat("✓ Results saved to output/regional_robustness_results.RData\n")
