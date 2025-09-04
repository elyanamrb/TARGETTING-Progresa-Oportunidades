# ============================================================================
# 01_load_and_clean_data.R
# ============================================================================
# Purpose: Load and clean all raw data files for Progresa-Oportunidades RDD analysis
# Author: [Your Name]
# Date Created: 2024
# Last Modified: 2024
#
# Inputs:
#   - iml1995.dta: Marginalization index 1995
#   - iml2000.dta: Marginalization index 2000  
#   - iml2005.dta: Marginalization index 2005
#   - entry_loc_300710.dta: Program entry years (source 1)
#   - entry_locality_enca.dta: Program entry years (source 2)
#   - regiones41.dta: Regional classifications
#   - Base_marginacion_localidad_90-10.csv: Panel data 1990-2010
#
# Outputs:
#   - Clean RDS files for each dataset
#   - Data quality report
# ============================================================================

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(haven)      # For reading Stata files
library(here)       # For project-relative paths
library(janitor)    # For data cleaning utilities
library(skimr)      # For data summaries
library(assertr)    # For data validation

# Set up paths (adjust base path as needed)
base_path <- "~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS"
setwd(base_path)

# Create directory structure if it doesn't exist
dir.create("issue-1.1-data-preparation", showWarnings = FALSE)
dir.create("issue-1.1-data-preparation/data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("issue-1.1-data-preparation/data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("issue-1.1-data-preparation/output/reports", recursive = TRUE, showWarnings = FALSE)
dir.create("issue-1.1-data-preparation/output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("issue-1.1-data-preparation/code", recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# SECTION 1: LOAD RAW DATA
# ============================================================================

message("Loading raw data files...")

# 1.1 Load IML datasets
iml_95_raw <- read_dta("iml1995.dta")
iml_00_raw <- read_dta("iml2000.dta")
iml_05_raw <- read_dta("iml2005.dta")

# 1.2 Load program entry datasets
entry_30_raw <- read_dta("entry_loc_300710.dta")
entry_enca_raw <- read_dta("entry_locality_enca.dta")

# 1.3 Load regional classification
regiones41_raw <- read_dta("regiones41.dta")

# 1.4 Load panel data (with encoding handling)
iml_panel_raw <- read.csv("Base_marginacion_localidad_90-10.csv", 
                          fileEncoding = "latin1",
                          stringsAsFactors = FALSE)

# Convert column names to lowercase for consistency
names(iml_panel_raw) <- tolower(names(iml_panel_raw))

# ============================================================================
# SECTION 2: CLEAN AND STANDARDIZE IML DATASETS
# ============================================================================

message("Cleaning IML datasets...")

# 2.1 Clean IML 1995
iml_95_clean <- iml_95_raw %>%
  # Standardize variable names
  rename(
    cve_edo = cve_ent,              # State code
    pobtot_95 = pobtot95,           # Total population 1995
    analf_95 = analf95,             # Illiteracy rate 1995
    sinag_95 = sinag95,             # Without water 1995
    sindr_95 = sindr95,             # Without drainage 1995
    sinel_95 = sinel95,             # Without electricity 1995
    pisotier_95 = pisotier90,       # Dirt floor 1995
    iml_95 = indice95               # Marginalization index 1995
  ) %>%
  # Select relevant variables
  select(cve_edo, cve_mun, cve_loc, iml_95, pobtot_95, 
         analf_95, sinag_95, sinel_95, pisotier_95) %>%
  # Convert to numeric
  mutate(
    iml_95 = as.numeric(iml_95),
    pobtot_95 = as.numeric(pobtot_95)
  ) %>%
  # Create unique locality ID
  mutate(
    loc_id = paste0(
      sprintf("%02d", cve_edo),
      sprintf("%03d", cve_mun),
      sprintf("%04d", cve_loc)
    )
  )

# 2.2 Clean IML 2000
iml_00_clean <- iml_00_raw %>%
  rename(
    cve_edo = cve_ent,
    pobtot_00 = pobtot00,
    analf_00 = anal_iml,
    sinag_00 = agua_iml,
    sinel_00 = elec_iml,
    pisotier_00 = tier_iml,
    iml_00 = imlr2000
  ) %>%
  select(cve_edo, cve_mun, cve_loc, iml_00, pobtot_00,
         analf_00, sinag_00, sinel_00, pisotier_00) %>%
  mutate(
    iml_00 = as.numeric(iml_00),
    pobtot_00 = as.numeric(pobtot_00),
    loc_id = paste0(
      sprintf("%02d", cve_edo),
      sprintf("%03d", cve_mun),
      sprintf("%04d", cve_loc)
    )
  )

# 2.3 Clean IML 2005
iml_05_clean <- iml_05_raw %>%
  rename(
    cve_edo = cve_ent,
    pobtot_05 = pobtot05,
    analf_05 = anal05,
    sinag_05 = sagua05,
    sindr_05 = sdye05,
    sinel_05 = see05,
    pisotier_05 = pisotie05,
    sinrefri_05 = srefri05,
    iml_05 = im_05
  ) %>%
  select(cve_edo, cve_mun, cve_loc, iml_05, pobtot_05,
         analf_05, sinag_05, sinel_05, pisotier_05, sinrefri_05) %>%
  mutate(
    iml_05 = as.numeric(iml_05),
    pobtot_05 = as.numeric(pobtot_05),
    loc_id = paste0(
      sprintf("%02d", cve_edo),
      sprintf("%03d", cve_mun),
      sprintf("%04d", cve_loc)
    )
  )

# ============================================================================
# SECTION 3: CLEAN PROGRAM ENTRY DATA
# ============================================================================

message("Cleaning program entry data...")

# 3.1 Clean first entry dataset
entry_30_clean <- entry_30_raw %>%
  rename(entry_year_30 = entry_year) %>%
  mutate(
    # Standardize municipality codes to 3 digits
    cve_mun_std = sprintf("%03d", cve_mun),
    # Create combined state-municipality code
    cve_mun_full = paste0(cve_edo, cve_mun_std),
    # Create unique locality ID
    loc_id = paste0(
      sprintf("%02d", cve_edo),
      sprintf("%03d", cve_mun),
      sprintf("%04d", cve_loc)
    )
  )

# 3.2 Clean second entry dataset
entry_enca_clean <- entry_enca_raw %>%
  rename(entry_year_enca = entry_year_enca) %>%
  mutate(
    loc_id = paste0(
      sprintf("%02d", cve_edo),
      sprintf("%03d", cve_mun),
      sprintf("%04d", cve_loc)
    )
  )

# ============================================================================
# SECTION 4: CLEAN REGIONAL CLASSIFICATION
# ============================================================================

message("Cleaning regional classification...")

regiones41_clean <- regiones41_raw %>%
  select(cve_edo, cve_mun, region41) %>%
  mutate(
    # Standardize municipality codes
    cve_mun_std = sprintf("%03d", cve_mun),
    cve_mun_full = paste0(cve_edo, cve_mun_std)
  )

# ============================================================================
# SECTION 5: CLEAN PANEL DATA
# ============================================================================

message("Cleaning panel data...")

iml_panel_clean <- iml_panel_raw %>%
  rename(
    cve_edo = cve_ent,
    year = año,  # Note: May need adjustment based on actual column name
    mge_id = cve_loc
  ) %>%
  # Remove any NA years
  filter(!is.na(year)) %>%
  # Ensure numeric types
  mutate(across(starts_with("iml"), as.numeric))

# ============================================================================
# SECTION 6: DATA VALIDATION
# ============================================================================

message("Running data validation checks...")

# Validation function
validate_dataset <- function(data, dataset_name) {
  
  message(paste("Validating", dataset_name, "..."))
  
  # Check for duplicates
  n_dups <- data %>% 
    group_by(loc_id) %>% 
    filter(n() > 1) %>% 
    nrow()
  
  if(n_dups > 0) {
    warning(paste(dataset_name, "has", n_dups, "duplicate localities"))
  }
  
  # Check for missing values in key variables
  missing_summary <- data %>%
    select(any_of(c("loc_id", "cve_edo", "cve_mun", "cve_loc"))) %>%
    summarise_all(~sum(is.na(.)))
  
  print(paste(dataset_name, "missing values:"))
  print(missing_summary)
  
  return(invisible(NULL))
}

# Run validations
validate_dataset(iml_95_clean, "IML 1995")
validate_dataset(iml_00_clean, "IML 2000")
validate_dataset(iml_05_clean, "IML 2005")
validate_dataset(entry_30_clean, "Entry 30")
validate_dataset(entry_enca_clean, "Entry ENCA")

# ============================================================================
# SECTION 7: SAVE CLEAN DATASETS
# ============================================================================

message("Saving clean datasets...")

# Save as RDS for preserving data types
saveRDS(iml_95_clean, "issue-1.1-data-preparation/data/processed/iml_95_clean.rds")
saveRDS(iml_00_clean, "issue-1.1-data-preparation/data/processed/iml_00_clean.rds")
saveRDS(iml_05_clean, "issue-1.1-data-preparation/data/processed/iml_05_clean.rds")
saveRDS(entry_30_clean, "issue-1.1-data-preparation/data/processed/entry_30_clean.rds")
saveRDS(entry_enca_clean, "issue-1.1-data-preparation/data/processed/entry_enca_clean.rds")
saveRDS(regiones41_clean, "issue-1.1-data-preparation/data/processed/regiones41_clean.rds")
saveRDS(iml_panel_clean, "issue-1.1-data-preparation/data/processed/iml_panel_clean.rds")

# ============================================================================
# SECTION 8: GENERATE SUMMARY REPORT
# ============================================================================

message("Generating data quality report...")

# Create summary statistics for each dataset
summary_list <- list(
  "IML 1995" = skim(iml_95_clean),
  "IML 2000" = skim(iml_00_clean),
  "IML 2005" = skim(iml_05_clean),
  "Entry 30" = skim(entry_30_clean),
  "Entry ENCA" = skim(entry_enca_clean),
  "Regions" = skim(regiones41_clean)
)

# Save summary report
sink("issue-1.1-data-preparation/output/reports/data_quality_report.txt")
cat("DATA QUALITY REPORT\n")
cat("===================\n")
cat("Generated:", Sys.Date(), "\n\n")

for(name in names(summary_list)) {
  cat("\n", name, "\n")
  cat(rep("-", nchar(name)), "\n")
  print(summary_list[[name]])
}
sink()

message("Step 1: Data loading and cleaning complete!")
message("Clean datasets saved to: issue-1.1-data-preparation/data/processed/")
message("Next step: Run 02_merge_datasets.R")
