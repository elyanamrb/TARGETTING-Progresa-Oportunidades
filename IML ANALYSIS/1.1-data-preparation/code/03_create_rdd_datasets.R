
# ============================================================================
# 03_create_rdd_datasets.R
# ============================================================================
# Purpose: Create specialized datasets for RDD analysis with optimal cutoffs
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# Inputs:
#   - merged_entry_iml_rg.csv from 02_merge_datasets.R
#
# Outputs:
#   - regional_cutoffs.csv: Optimal cutoffs by region
#   - rdd_analysis_dataset.csv: Dataset with normalized running variable
#   - Summary statistics tables
# ============================================================================

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(rdrobust)  # For RDD optimal bandwidth selection
library(rddensity) # For manipulation testing
library(broom)     # For tidy model outputs

# Set working directory
base_path <- "~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS"
setwd(base_path)

# ============================================================================
# SECTION 1: LOAD MERGED DATA
# ============================================================================

message("Loading merged dataset...")

df <- read_csv("issue-1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

# ============================================================================
# SECTION 2: FIND OPTIMAL CUTOFFS BY REGION
# ============================================================================

message("Finding optimal cutoffs by region...")

# Function to find optimal cutoff for a region
find_optimal_cutoff <- function(data, region_id) {
  
  # Filter to region
  region_data <- data %>%
    filter(region41 == region_id, !is.na(iml_95))
  
  if(nrow(region_data) < 50) {
    return(data.frame(
      region = region_id,
      optimal_cutoff = NA,
      n_obs = nrow(region_data),
      rd_estimate = NA,
      se = NA,
      pvalue = NA,
      bandwidth = NA
    ))
  }
  
  # Try different potential cutoffs (based on quantiles)
  potential_cutoffs <- quantile(region_data$iml_95, 
                                probs = seq(0.1, 0.6, by = 0.05),
                                na.rm = TRUE)
  
  results <- map_df(potential_cutoffs, function(cutoff) {
    
    # Run RDD at this cutoff
    tryCatch({
      rd_result <- rdrobust(
        y = region_data$post_2002,
        x = region_data$iml_95,
        c = cutoff,
        kernel = "triangular",
        bwselect = "mserd"
      )
      
      return(data.frame(
        cutoff = cutoff,
        estimate = rd_result$coef[1],
        se = rd_result$se[1],
        pvalue = rd_result$pv[1],
        bandwidth = rd_result$bws[1],
        n_left = sum(region_data$iml_95 < cutoff),
        n_right = sum(region_data$iml_95 >= cutoff)
      ))
    }, error = function(e) {
      return(NULL)
    })
  })
  
  if(nrow(results) == 0) {
    return(data.frame(
      region = region_id,
      optimal_cutoff = NA,
      n_obs = nrow(region_data),
      rd_estimate = NA,
      se = NA,
      pvalue = NA,
      bandwidth = NA
    ))
  }
  
  # Select cutoff with strongest (most negative) effect
  optimal <- results %>%
    filter(pvalue < 0.10) %>%  # Only consider significant results
    arrange(estimate) %>%
    slice(1)
  
  if(nrow(optimal) == 0) {
    optimal <- results %>% arrange(estimate) %>% slice(1)
  }
  
  return(data.frame(
    region = region_id,
    optimal_cutoff = optimal$cutoff,
    n_obs = nrow(region_data),
    rd_estimate = optimal$estimate,
    se = optimal$se,
    pvalue = optimal$pvalue,
    bandwidth = optimal$bandwidth
  ))
}

# Find optimal cutoffs for all regions
regions <- unique(df$region41[!is.na(df$region41)])
regional_cutoffs <- map_df(regions, ~find_optimal_cutoff(df, .x))

# Add region labels
regional_cutoffs <- regional_cutoffs %>%
  arrange(region) %>%
  mutate(
    significance = case_when(
      pvalue < 0.01 ~ "***",
      pvalue < 0.05 ~ "**",
      pvalue < 0.10 ~ "*",
      TRUE ~ ""
    )
  )

# Save regional cutoffs
write_csv(regional_cutoffs, 
          "issue-1.1-data-preparation/data/processed/regional_cutoffs.csv")

message(paste("Optimal cutoffs found for", 
              sum(!is.na(regional_cutoffs$optimal_cutoff)), 
              "out of", 
              nrow(regional_cutoffs), 
              "regions"))

# ============================================================================
# SECTION 3: CREATE NORMALIZED RDD DATASET
# ============================================================================

message("Creating normalized RDD dataset...")

# Merge cutoffs with main data
df_with_cutoffs <- df %>%
  left_join(
    regional_cutoffs %>% select(region = region, optimal_cutoff),
    by = c("region41" = "region")
  ) %>%
  mutate(
    # Create normalized running variable (distance from regional cutoff)
    iml_95_normalized = iml_95 - optimal_cutoff,
    
    # Create treatment indicator based on cutoff
    above_cutoff = ifelse(iml_95_normalized > 0, 1, 0),
    
    # Distance to cutoff (absolute value)
    dist_to_cutoff = abs(iml_95_normalized)
  )

# ============================================================================
# SECTION 4: CREATE ANALYSIS SAMPLES
# ============================================================================

message("Creating analysis samples...")

# Main analysis sample (localities with valid cutoffs)
rdd_main <- df_with_cutoffs %>%
  filter(!is.na(optimal_cutoff))

# Bandwidth-restricted samples
rdd_narrow <- rdd_main %>%
  filter(abs(iml_95_normalized) <= 0.5)  # Narrow bandwidth

rdd_wide <- rdd_main %>%
  filter(abs(iml_95_normalized) <= 1.5)  # Wide bandwidth

# Save RDD datasets
write_csv(rdd_main, 
          "issue-1.1-data-preparation/data/processed/rdd_analysis_dataset.csv")
write_csv(rdd_narrow,
          "issue-1.1-data-preparation/data/processed/rdd_narrow_bandwidth.csv")
write_csv(rdd_wide,
          "issue-1.1-data-preparation/data/processed/rdd_wide_bandwidth.csv")

# ============================================================================
# SECTION 5: GENERATE SUMMARY STATISTICS
# ============================================================================

message("Generating summary statistics...")

# Function to create summary statistics table
create_summary_table <- function(data, sample_name) {
  
  # Variables to summarize
  vars_continuous <- c("iml_95", "iml_00", "iml_05", "pobtot_95", 
                       "analf_95", "sinag_95", "sinel_95", "pisotier_95")
  
  # Select available variables and convert to numeric
  data_for_summary <- data %>%
    select(all_of(vars_continuous[vars_continuous %in% names(data)])) %>%
    mutate(across(everything(), as.numeric))  # Convert all to numeric
  
  # Summary for continuous variables
  summary_stats <- data_for_summary %>%
    pivot_longer(everything(), names_to = "variable") %>%
    group_by(variable) %>%
    summarise(
      n = sum(!is.na(value)),
      mean = round(mean(value, na.rm = TRUE), 3),
      sd = round(sd(value, na.rm = TRUE), 3),
      min = round(min(value, na.rm = TRUE), 3),
      p25 = round(quantile(value, 0.25, na.rm = TRUE), 3),
      median = round(median(value, na.rm = TRUE), 3),
      p75 = round(quantile(value, 0.75, na.rm = TRUE), 3),
      max = round(max(value, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    mutate(sample = sample_name)
  
  return(summary_stats)
}

# Create summary tables for different samples
summary_full <- create_summary_table(df, "Full Sample")
summary_rdd <- create_summary_table(rdd_main, "RDD Sample")
summary_pre2002 <- create_summary_table(
  rdd_main %>% filter(post_2002 == 0), 
  "Pre-2002 Entry"
)
summary_post2002 <- create_summary_table(
  rdd_main %>% filter(post_2002 == 1), 
  "Post-2002 Entry"
)

# Combine all summaries
all_summaries <- bind_rows(
  summary_full,
  summary_rdd,
  summary_pre2002,
  summary_post2002
)

# Save summary statistics
write_csv(all_summaries,
          "issue-1.1-data-preparation/output/tables/summary_statistics.csv")

# ============================================================================
# SECTION 6: BALANCE TESTS
# ============================================================================

message("Performing balance tests...")

# Ensure all variables are numeric before balance tests
rdd_main <- rdd_main %>%
  mutate(across(c(pobtot_95, analf_95, sinag_95, sinel_95, pisotier_95), as.numeric))

# Test balance at cutoff for pre-treatment variables
balance_tests <- rdd_main %>%
  filter(abs(iml_95_normalized) <= 0.5) %>%  # Use narrow bandwidth
  select(above_cutoff, pobtot_95, analf_95, sinag_95, sinel_95, pisotier_95) %>%
  pivot_longer(-above_cutoff, names_to = "variable") %>%
  group_by(variable) %>%
  nest() %>%
  mutate(
    test = map(data, ~t.test(value ~ above_cutoff, data = .x)),
    tidy_test = map(test, tidy)
  ) %>%
  unnest(tidy_test) %>%
  select(variable, estimate1, estimate2, estimate, p.value) %>%
  rename(
    mean_below = estimate1,
    mean_above = estimate2,
    difference = estimate
  )

# Save balance tests
write_csv(balance_tests,
          "issue-1.1-data-preparation/output/tables/balance_tests.csv")

# ============================================================================
# SECTION 7: MANIPULATION TESTING
# ============================================================================

message("Testing for manipulation at cutoffs...")

# Test for manipulation at cutoff (McCrary test)
manipulation_test <- tryCatch({
  rddensity(rdd_main$iml_95_normalized, c = 0)
}, error = function(e) {
  message("Manipulation test failed: ", e$message)
  return(NULL)
})

if(!is.null(manipulation_test)) {
  # Save test results
  sink("issue-1.1-data-preparation/output/reports/manipulation_test.txt")
  cat("MANIPULATION TEST RESULTS\n")
  cat("=========================\n")
  cat("Test statistic:", manipulation_test$test$t_jk, "\n")
  cat("P-value:", manipulation_test$test$p_jk, "\n")
  cat("\nInterpretation:\n")
  if(manipulation_test$test$p_jk > 0.05) {
    cat("No evidence of manipulation at the cutoff (p > 0.05)\n")
  } else {
    cat("Warning: Potential manipulation detected at cutoff (p < 0.05)\n")
  }
  sink()
}

# ============================================================================
# SECTION 8: ENTRY PATTERN ANALYSIS
# ============================================================================

message("Analyzing program entry patterns...")

# Create entry pattern summary
entry_patterns <- df %>%
  group_by(entry_year) %>%
  summarise(
    n_localities = n(),
    mean_iml_95 = mean(iml_95, na.rm = TRUE),
    sd_iml_95 = sd(iml_95, na.rm = TRUE),
    median_iml_95 = median(iml_95, na.rm = TRUE),
    min_iml_95 = min(iml_95, na.rm = TRUE),
    max_iml_95 = max(iml_95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_total = n_localities / sum(n_localities) * 100,
    cumulative_n = cumsum(n_localities),
    cumulative_pct = cumulative_n / sum(n_localities) * 100,
    program_phase = ifelse(entry_year <= 2002, "Progresa", "Oportunidades")
  )

write_csv(entry_patterns,
          "issue-1.1-data-preparation/output/tables/entry_patterns.csv")

# Regional entry patterns
regional_entry_patterns <- df %>%
  group_by(region41, entry_year) %>%
  summarise(
    n_localities = n(),
    mean_iml_95 = mean(iml_95, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(region41) %>%
  mutate(
    pct_regional = n_localities / sum(n_localities) * 100
  ) %>%
  ungroup()

write_csv(regional_entry_patterns,
          "issue-1.1-data-preparation/output/tables/regional_entry_patterns.csv")

# Identify peak entry years
peak_years <- entry_patterns %>%
  arrange(desc(n_localities)) %>%
  slice(1:5) %>%
  select(entry_year, n_localities, pct_total, program_phase)

message("Peak entry years:")
print(peak_years)

# ============================================================================
# SECTION 9: FINAL REPORT
# ============================================================================

message("Generating final data preparation report...")

sink("issue-1.1-data-preparation/output/reports/data_preparation_summary.txt")
cat("DATA PREPARATION SUMMARY\n")
cat("========================\n")
cat("Generated:", Sys.Date(), "\n\n")

cat("DATASET DIMENSIONS\n")
cat("------------------\n")
cat("Original merged dataset:", nrow(df), "localities\n")
cat("RDD analysis dataset:", nrow(rdd_main), "localities\n")
cat("Narrow bandwidth sample:", nrow(rdd_narrow), "localities\n")
cat("Wide bandwidth sample:", nrow(rdd_wide), "localities\n\n")

cat("REGIONAL COVERAGE\n")
cat("-----------------\n")
cat("Total regions:", n_distinct(df$region41, na.rm = TRUE), "\n")
cat("Regions with valid cutoffs:", sum(!is.na(regional_cutoffs$optimal_cutoff)), "\n")
cat("Average localities per region:", round(mean(regional_cutoffs$n_obs)), "\n\n")

cat("TREATMENT STATUS\n")
cat("----------------\n")
cat("Pre-2002 entries:", sum(rdd_main$post_2002 == 0), 
    "(", round(mean(rdd_main$post_2002 == 0) * 100, 1), "%)\n")
cat("Post-2002 entries:", sum(rdd_main$post_2002 == 1),
    "(", round(mean(rdd_main$post_2002 == 1) * 100, 1), "%)\n\n")

cat("CUTOFF STATISTICS\n")
cat("-----------------\n")
cat("Mean optimal cutoff:", round(mean(regional_cutoffs$optimal_cutoff, na.rm = TRUE), 3), "\n")
cat("SD of cutoffs:", round(sd(regional_cutoffs$optimal_cutoff, na.rm = TRUE), 3), "\n")
cat("Range:", round(min(regional_cutoffs$optimal_cutoff, na.rm = TRUE), 3), 
    "to", round(max(regional_cutoffs$optimal_cutoff, na.rm = TRUE), 3), "\n")

sink()

message("\nStep 3: RDD dataset creation complete!")
message("Files saved to: issue-1.1-data-preparation/data/processed/")
message("Reports saved to: issue-1.1-data-preparation/output/reports/")
message("Tables saved to: issue-1.1-data-preparation/output/tables/")
message("\nData preparation pipeline complete! Ready for RDD analysis.")
