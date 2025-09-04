# =============================================================================
# Script: 02_prepare_and_harmonize_data.R
# Purpose: Prepare and harmonize ENCASEH data across years for analysis
#          - Standardize variable coding
#          - Handle asset variable inconsistencies
#          - Create analysis-ready variables
#          - Handle missing values and outliers
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(haven)
library(data.table)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/01-data-foundation")

# Create output directories
dir.create("outputs/data_quality", showWarnings = FALSE)

cat("========================================\n")
cat("DATA PREPARATION AND HARMONIZATION\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD EXPLORATION RESULTS
# =============================================================================

cat("Loading exploration results...\n")
exploration_data <- readRDS("data/processed/exploration_results.rds")

# Extract individual years
encaseh_2000 <- exploration_data$data_years$y2000
encaseh_2001 <- exploration_data$data_years$y2001
encaseh_2002 <- exploration_data$data_years$y2002
encaseh_2004 <- exploration_data$data_years$y2004

# =============================================================================
# PART 2: STANDARDIZE ASSET VARIABLES
# =============================================================================

cat("\n========================================\n")
cat("STANDARDIZING ASSET VARIABLES\n")
cat("========================================\n")

# Define asset variables
asset_vars <- c("LICUADORA", "REFRI", "EST_GAS", "CAL_AGUA", "RADIO",
                "TOCADISC", "TELEVIS", "VIDEOCAS", "LAVADORA", "VENTILA",
                "AUTO_PROP", "CAMION_PRO", "TRACTOR")

# Function to standardize asset coding
# Based on your documentation: 1 = has asset, 2 = doesn't have asset
# Need to handle: 0, -1, 9, NA
standardize_assets <- function(data, year) {
  cat(sprintf("\nStandardizing assets for year %d...\n", year))
  
  # Find which asset variables are present
  present_assets <- intersect(asset_vars, names(data))
  
  for (asset in present_assets) {
    # Check unique values before standardization
    unique_vals <- unique(data[[asset]])
    cat(sprintf("  %s unique values: %s\n", asset, 
                paste(sort(unique_vals[!is.na(unique_vals)]), collapse = ", ")))
    
    # Standardize to binary: 1 = has, 0 = doesn't have
    data[[paste0(asset, "_std")]] <- case_when(
      data[[asset]] == 1 ~ 1,  # Has asset
      data[[asset]] == 2 ~ 0,  # Doesn't have asset
      data[[asset]] == 0 ~ 0,  # Treat 0 as doesn't have
      data[[asset]] == 9 ~ NA_real_,  # 9 likely missing/refused
      data[[asset]] == -1 ~ NA_real_, # -1 likely missing
      TRUE ~ NA_real_
    )
    
    # Report standardization results
    n_has <- sum(data[[paste0(asset, "_std")]] == 1, na.rm = TRUE)
    n_hasnt <- sum(data[[paste0(asset, "_std")]] == 0, na.rm = TRUE)
    n_missing <- sum(is.na(data[[paste0(asset, "_std")]]))
    
    cat(sprintf("    Standardized: %d have, %d don't have, %d missing\n",
                n_has, n_hasnt, n_missing))
  }
  
  # Calculate asset index (number of assets owned)
  asset_std_vars <- paste0(present_assets, "_std")
  data$asset_count <- rowSums(data[asset_std_vars], na.rm = TRUE)
  data$asset_index <- data$asset_count / length(present_assets)
  
  cat(sprintf("  Created asset index: mean = %.2f assets\n", 
              mean(data$asset_count, na.rm = TRUE)))
  
  return(data)
}

# Standardize each year
if (!is.null(encaseh_2000)) encaseh_2000 <- standardize_assets(encaseh_2000, 2000)
if (!is.null(encaseh_2001)) encaseh_2001 <- standardize_assets(encaseh_2001, 2001)
if (!is.null(encaseh_2002)) encaseh_2002 <- standardize_assets(encaseh_2002, 2002)
if (!is.null(encaseh_2004)) encaseh_2004 <- standardize_assets(encaseh_2004, 2004)

# =============================================================================
# PART 3: HARMONIZE TARGETING METHOD (TIPO_CAP)
# =============================================================================

cat("\n========================================\n")
cat("HARMONIZING TARGETING METHODS\n")
cat("========================================\n")

harmonize_targeting <- function(data, year) {
  cat(sprintf("\nProcessing targeting method for year %d...\n", year))
  
  if ("TIPO_CAP" %in% names(data)) {
    # Create readable targeting method variable
    data$targeting_method <- case_when(
      data$TIPO_CAP == 4 ~ "Ventanilla",
      data$TIPO_CAP == 5 ~ "Rural Census",
      data$TIPO_CAP == 6 ~ "Urban Census",
      TRUE ~ NA_character_
    )
    
    # Create binary indicators
    data$is_ventanilla <- ifelse(data$TIPO_CAP == 4, 1, 0)
    data$is_census <- ifelse(data$TIPO_CAP %in% c(5, 6), 1, 0)
    data$is_rural <- ifelse(data$TIPO_CAP == 5, 1, 0)
    data$is_urban <- ifelse(data$TIPO_CAP == 6, 1, 0)
    
    # Summary
    method_table <- table(data$targeting_method, useNA = "ifany")
    cat("  Targeting method distribution:\n")
    print(method_table)
    
  } else {
    cat(sprintf("  TIPO_CAP not found in year %d data\n", year))
    # Add NA columns for consistency
    data$targeting_method <- NA_character_
    data$is_ventanilla <- NA
    data$is_census <- NA
    data$is_rural <- NA
    data$is_urban <- NA
  }
  
  return(data)
}

# Harmonize each year
if (!is.null(encaseh_2000)) encaseh_2000 <- harmonize_targeting(encaseh_2000, 2000)
if (!is.null(encaseh_2001)) encaseh_2001 <- harmonize_targeting(encaseh_2001, 2001)
if (!is.null(encaseh_2002)) encaseh_2002 <- harmonize_targeting(encaseh_2002, 2002)
if (!is.null(encaseh_2004)) encaseh_2004 <- harmonize_targeting(encaseh_2004, 2004)

# =============================================================================
# PART 4: STANDARDIZE OUTCOME VARIABLE (RESULTADO)
# =============================================================================

cat("\n========================================\n")
cat("STANDARDIZING OUTCOME VARIABLE\n")
cat("========================================\n")

standardize_outcome <- function(data, year) {
  cat(sprintf("\nStandardizing RESULTADO for year %d...\n", year))
  
  if ("RESULTADO" %in% names(data)) {
    # Check unique values
    unique_vals <- unique(data$RESULTADO)
    cat(sprintf("  Unique RESULTADO values: %s\n", 
                paste(sort(unique_vals[!is.na(unique_vals)]), collapse = ", ")))
    
    # Create binary accepted variable
    # Based on analysis: RESULTADO = 1 likely means accepted
    data$accepted <- ifelse(data$RESULTADO == 1, 1, 0)
    
    # Summary
    accept_rate <- mean(data$accepted, na.rm = TRUE) * 100
    cat(sprintf("  Acceptance rate: %.1f%%\n", accept_rate))
    
  } else {
    cat(sprintf("  RESULTADO not found in year %d\n", year))
    data$accepted <- NA
  }
  
  return(data)
}

# Standardize each year
if (!is.null(encaseh_2000)) encaseh_2000 <- standardize_outcome(encaseh_2000, 2000)
if (!is.null(encaseh_2001)) encaseh_2001 <- standardize_outcome(encaseh_2001, 2001)
if (!is.null(encaseh_2002)) encaseh_2002 <- standardize_outcome(encaseh_2002, 2002)
if (!is.null(encaseh_2004)) encaseh_2004 <- standardize_outcome(encaseh_2004, 2004)

# =============================================================================
# PART 5: HANDLE POVERTY SCORE (PUN_INC)
# =============================================================================

cat("\n========================================\n")
cat("PROCESSING POVERTY SCORES\n")
cat("========================================\n")

process_poverty_score <- function(data, year) {
  cat(sprintf("\nProcessing PUN_INC for year %d...\n", year))
  
  if ("PUN_INC" %in% names(data)) {
    # Identify scoring system
    data$scoring_system <- ifelse(data$PUN_INC > 100, "System_2", "System_1")
    
    # Create cleaned version (for analyses that need it)
    data$pun_inc_clean <- ifelse(data$PUN_INC > 100, NA, data$PUN_INC)
    
    # Create poverty quintiles and deciles (using clean scores)
    data$poverty_quintile <- ntile(data$pun_inc_clean, 5)
    data$poverty_decile <- ntile(data$pun_inc_clean, 10)
    
    # Summary
    system_table <- table(data$scoring_system)
    cat("  Scoring system distribution:\n")
    print(system_table)
    
    cat(sprintf("  Clean PUN_INC: mean = %.3f, median = %.3f\n",
                mean(data$pun_inc_clean, na.rm = TRUE),
                median(data$pun_inc_clean, na.rm = TRUE)))
  }
  
  return(data)
}

# Process each year
if (!is.null(encaseh_2000)) encaseh_2000 <- process_poverty_score(encaseh_2000, 2000)
if (!is.null(encaseh_2001)) encaseh_2001 <- process_poverty_score(encaseh_2001, 2001)
if (!is.null(encaseh_2002)) encaseh_2002 <- process_poverty_score(encaseh_2002, 2002)
if (!is.null(encaseh_2004)) encaseh_2004 <- process_poverty_score(encaseh_2004, 2004)

# =============================================================================
# PART 6: ADD YEAR IDENTIFIER
# =============================================================================

cat("\n========================================\n")
cat("ADDING YEAR IDENTIFIERS\n")
cat("========================================\n")

if (!is.null(encaseh_2000)) encaseh_2000$year <- 2000
if (!is.null(encaseh_2001)) encaseh_2001$year <- 2001
if (!is.null(encaseh_2002)) encaseh_2002$year <- 2002
if (!is.null(encaseh_2004)) encaseh_2004$year <- 2004

# =============================================================================
# PART 7: CREATE DATA QUALITY REPORT
# =============================================================================

cat("\n========================================\n")
cat("GENERATING DATA QUALITY REPORT\n")
cat("========================================\n")

create_quality_report <- function(data, year) {
  if (is.null(data)) return(NULL)
  
  report <- data.frame(
    year = year,
    n_obs = nrow(data),
    n_accepted = sum(data$accepted == 1, na.rm = TRUE),
    accept_rate = mean(data$accepted, na.rm = TRUE) * 100,
    n_ventanilla = sum(data$is_ventanilla == 1, na.rm = TRUE),
    n_census = sum(data$is_census == 1, na.rm = TRUE),
    mean_assets = mean(data$asset_count, na.rm = TRUE),
    pct_system2 = mean(data$scoring_system == "System_2", na.rm = TRUE) * 100,
    n_missing_pun_inc = sum(is.na(data$PUN_INC)),
    n_missing_resultado = sum(is.na(data$RESULTADO))
  )
  
  return(report)
}

# Create quality reports
quality_2000 <- create_quality_report(encaseh_2000, 2000)
quality_2001 <- create_quality_report(encaseh_2001, 2001)
quality_2002 <- create_quality_report(encaseh_2002, 2002)
quality_2004 <- create_quality_report(encaseh_2004, 2004)

quality_summary <- bind_rows(quality_2000, quality_2001, quality_2002, quality_2004)

# Save quality report
write.csv(quality_summary, "outputs/data_quality/quality_summary.csv", row.names = FALSE)

cat("\nData Quality Summary:\n")
print(quality_summary)

# =============================================================================
# PART 8: SAVE HARMONIZED DATA
# =============================================================================

cat("\n========================================\n")
cat("SAVING HARMONIZED DATA\n")
cat("========================================\n")

# Save individual years
if (!is.null(encaseh_2000)) {
  saveRDS(encaseh_2000, "data/processed/encaseh_2000_harmonized.rds")
  cat("✓ Saved 2000 data\n")
}

if (!is.null(encaseh_2001)) {
  saveRDS(encaseh_2001, "data/processed/encaseh_2001_harmonized.rds")
  cat("✓ Saved 2001 data\n")
}

if (!is.null(encaseh_2002)) {
  saveRDS(encaseh_2002, "data/processed/encaseh_2002_harmonized.rds")
  cat("✓ Saved 2002 data\n")
}

if (!is.null(encaseh_2004)) {
  saveRDS(encaseh_2004, "data/processed/encaseh_2004_harmonized.rds")
  cat("✓ Saved 2004 data\n")
}

# Save combined dataset for analyses that need it
all_years <- bind_rows(
  encaseh_2000,
  encaseh_2001,
  encaseh_2002,
  encaseh_2004
)

saveRDS(all_years, "data/processed/encaseh_all_years_harmonized.rds")
cat("✓ Saved combined dataset\n")

# =============================================================================
# SUMMARY
# =============================================================================

cat("\n========================================\n")
cat("HARMONIZATION COMPLETE\n")
cat("========================================\n")
cat("\nKey harmonization steps completed:\n")
cat("1. ✓ Asset variables standardized to binary (1=has, 0=doesn't have)\n")
cat("2. ✓ Targeting methods harmonized (Ventanilla/Rural Census/Urban Census)\n")
cat("3. ✓ Outcome variable standardized (accepted = 1/0)\n")
cat("4. ✓ Poverty scores processed (dual system identified)\n")
cat("5. ✓ Year identifiers added\n")
cat("6. ✓ Data quality report generated\n")
cat("\nHarmonized data saved to 'data/processed/'\n")
