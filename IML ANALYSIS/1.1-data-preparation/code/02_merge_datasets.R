# ============================================================================
# 02_merge_datasets.R
# ============================================================================
# Purpose: Merge clean datasets to create analytical files for RDD analysis
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# Inputs:
#   - Clean RDS files from 01_load_and_clean_data.R
#
# Outputs:
#   - merged_entry_iml_rg.csv: Main cross-sectional dataset
#   - iml_panel_90_10.csv: Panel dataset with entry timing
#   - Merge diagnostics report
# ============================================================================

# Clear environment
rm(list = ls())

# Load required packages
library(tidyverse)
library(janitor)
library(assertr)

# Set working directory
base_path <- "~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS"
setwd(base_path)

# ============================================================================
# SECTION 1: LOAD CLEAN DATASETS
# ============================================================================

message("Loading clean datasets...")

iml_95 <- readRDS("issue-1.1-data-preparation/data/processed/iml_95_clean.rds")
iml_00 <- readRDS("issue-1.1-data-preparation/data/processed/iml_00_clean.rds")
iml_05 <- readRDS("issue-1.1-data-preparation/data/processed/iml_05_clean.rds")
entry_30 <- readRDS("issue-1.1-data-preparation/data/processed/entry_30_clean.rds")
entry_enca <- readRDS("issue-1.1-data-preparation/data/processed/entry_enca_clean.rds")
regiones41 <- readRDS("issue-1.1-data-preparation/data/processed/regiones41_clean.rds")
iml_panel <- readRDS("issue-1.1-data-preparation/data/processed/iml_panel_clean.rds")

# ============================================================================
# SECTION 2: MERGE ENTRY DATASETS
# ============================================================================

message("Merging entry datasets...")

# First, merge the two entry sources
# Priority: entry_30 > entry_enca
entry_combined <- entry_30 %>%
  select(cve_edo, cve_mun, cve_loc, loc_id, entry_year_30, cve_mun_full) %>%
  full_join(
    entry_enca %>% select(cve_edo, cve_mun, cve_loc, loc_id, entry_year_enca),
    by = c("cve_edo", "cve_mun", "cve_loc", "loc_id")
  ) %>%
  mutate(
    # Use entry_30 if available, otherwise entry_enca
    entry_year = coalesce(entry_year_30, entry_year_enca),
    # Create source indicator
    entry_source = case_when(
      !is.na(entry_year_30) & !is.na(entry_year_enca) ~ "both",
      !is.na(entry_year_30) ~ "entry_30",
      !is.na(entry_year_enca) ~ "entry_enca",
      TRUE ~ "none"
    )
  )

# Diagnostic: Check entry year distribution
entry_year_summary <- entry_combined %>%
  group_by(entry_year) %>%
  summarise(
    n_localities = n(),
    pct_localities = n() / nrow(entry_combined) * 100
  ) %>%
  arrange(entry_year)

message(paste("Entry year range:", 
              min(entry_combined$entry_year, na.rm = TRUE), 
              "to", 
              max(entry_combined$entry_year, na.rm = TRUE)))

# ============================================================================
# SECTION 3: CREATE MAIN CROSS-SECTIONAL DATASET
# ============================================================================

message("Creating main cross-sectional dataset...")

# Start with entry data and progressively merge IML scores and regions
merged_cross_section <- entry_combined %>%
  # Add IML 1995
  left_join(
    iml_95 %>% select(loc_id, iml_95, pobtot_95, analf_95, 
                      sinag_95, sinel_95, pisotier_95),
    by = "loc_id"
  ) %>%
  # Add IML 2000
  left_join(
    iml_00 %>% select(loc_id, iml_00, pobtot_00, analf_00,
                      sinag_00, sinel_00, pisotier_00),
    by = "loc_id"
  ) %>%
  # Add IML 2005
  left_join(
    iml_05 %>% select(loc_id, iml_05, pobtot_05, analf_05,
                      sinag_05, sinel_05, pisotier_05, sinrefri_05),
    by = "loc_id"
  ) %>%
  # Add regional classification
  left_join(
    regiones41 %>% select(cve_edo, cve_mun, region41),
    by = c("cve_edo", "cve_mun")
  )

# ============================================================================
# SECTION 4: CREATE RDD VARIABLES
# ============================================================================

message("Creating RDD analysis variables...")

merged_cross_section <- merged_cross_section %>%
  mutate(
    # Create post-2002 dummy (main treatment indicator)
    post_2002 = ifelse(entry_year > 2002, 1, 0),
    
    # Create post-2003 dummy (alternative specification)
    post_2003 = ifelse(entry_year > 2003, 1, 0),
    
    # Create entry period categories
    entry_period = case_when(
      entry_year <= 1998 ~ "1997-1998 (Progresa Phase I)",
      entry_year <= 2000 ~ "1999-2000 (Progresa Phase II)",
      entry_year <= 2002 ~ "2001-2002 (Transition)",
      entry_year <= 2004 ~ "2003-2004 (Oportunidades Early)",
      entry_year <= 2006 ~ "2005-2006 (Oportunidades Expansion)",
      entry_year > 2006 ~ "2007+ (Late Entry)",
      is.na(entry_year) ~ "Never Entered"
    ),
    
    # Create program phase indicator
    program_phase = case_when(
      entry_year <= 2002 ~ "Progresa",
      entry_year > 2002 ~ "Oportunidades",
      is.na(entry_year) ~ "Not Enrolled"
    )
  )

# ============================================================================
# SECTION 5: FILTER AND CLEAN
# ============================================================================

message("Applying filters and final cleaning...")

# Filter to valid entry years and non-missing IML 1995
merged_final <- merged_cross_section %>%
  filter(
    !is.na(iml_95),                    # Must have 1995 marginalization score
    !is.na(entry_year),                 # Must have entry year
    entry_year >= 1995,                 # Entry after 1995
    entry_year <= 2010                  # Entry before 2010
  ) %>%
  # Ensure all percentage and IML variables are numeric
  mutate(
    across(c(iml_95, iml_00, iml_05,
             pobtot_95, pobtot_00, pobtot_05,
             analf_95, analf_00, analf_05,
             sinag_95, sinag_00, sinag_05,
             sinel_95, sinel_00, sinel_05,
             pisotier_95, pisotier_00, pisotier_05,
             sinrefri_05), 
           as.numeric)
  ) %>%
  # Sort by state, municipality, locality
  arrange(cve_edo, cve_mun, cve_loc)

# Report filtering results
n_before <- nrow(merged_cross_section)
n_after <- nrow(merged_final)
message(paste("Localities before filtering:", n_before))
message(paste("Localities after filtering:", n_after))
message(paste("Localities dropped:", n_before - n_after, 
              "(", round((n_before - n_after)/n_before * 100, 1), "%)", sep = ""))

# ============================================================================
# SECTION 6: CREATE PANEL DATASET
# ============================================================================

message("Creating panel dataset...")

# The panel data has a different coding structure:
# - cve_mun contains state+municipality (e.g., 1001 = state 1, mun 001)
# - mge_id is 8 digits without leading zeros

# First, extract state and municipality from panel's combined cve_mun
iml_panel <- iml_panel %>%
  mutate(
    # Extract state from combined cve_mun (first 1-2 digits)
    cve_edo_panel = floor(cve_mun / 1000),
    # Extract municipality (last 3 digits)
    cve_mun_panel = cve_mun %% 1000,
    # Standardize mge_id to 9 digits with leading zero if needed
    mge_id_std = case_when(
      nchar(as.character(mge_id)) == 8 ~ paste0("0", as.character(mge_id)),
      TRUE ~ as.character(mge_id)
    )
  )

# Prepare entry data with standardized ID
entry_for_panel <- entry_combined %>%
  mutate(
    # Ensure mge_id (loc_id) is 9 digits
    mge_id_std = loc_id
  ) %>%
  select(mge_id_std, entry_year)

# Merge using standardized IDs
iml_panel_final <- iml_panel %>%
  left_join(
    entry_for_panel,
    by = "mge_id_std"
  ) %>%
  # Add regional classification using extracted state and municipality
  left_join(
    regiones41 %>% select(cve_edo, cve_mun, region41),
    by = c("cve_edo_panel" = "cve_edo", "cve_mun_panel" = "cve_mun")
  ) %>%
  # Create treatment indicators
  mutate(
    post_2002 = ifelse(entry_year > 2002, 1, 0),
    post_2003 = ifelse(entry_year > 2003, 1, 0),
    treated = ifelse(!is.na(entry_year) & year >= entry_year, 1, 0)
  ) %>%
  # Clean up - keep original mge_id and remove temporary columns
  select(-cve_edo_panel, -cve_mun_panel, -mge_id_std) %>%
  # Sort
  arrange(mge_id, year)

# ============================================================================
# SECTION 7: SAVE FINAL DATASETS
# ============================================================================

message("Saving final analytical datasets...")

# Save cross-sectional dataset
write_csv(merged_final, 
          "issue-1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

# Save panel dataset
write_csv(iml_panel_final,
          "issue-1.1-data-preparation/data/processed/iml_panel_90_10.csv")

# Also save as RDS for preserving data types
saveRDS(merged_final,
        "issue-1.1-data-preparation/data/processed/merged_entry_iml_rg.rds")
saveRDS(iml_panel_final,
        "issue-1.1-data-preparation/data/processed/iml_panel_90_10.rds")

# ============================================================================
# SECTION 8: GENERATE MERGE DIAGNOSTICS
# ============================================================================

message("Generating merge diagnostics...")

# Create diagnostics report
sink("issue-1.1-data-preparation/output/reports/merge_diagnostics.txt")
cat("MERGE DIAGNOSTICS REPORT\n")
cat("========================\n")
cat("Generated:", Sys.Date(), "\n\n")

cat("\n1. ENTRY DATA MERGE\n")
cat("-------------------\n")
cat("Localities with entry_30 only:", 
    sum(entry_combined$entry_source == "entry_30", na.rm = TRUE), "\n")
cat("Localities with entry_enca only:", 
    sum(entry_combined$entry_source == "entry_enca", na.rm = TRUE), "\n")
cat("Localities with both sources:", 
    sum(entry_combined$entry_source == "both", na.rm = TRUE), "\n")
cat("Localities with no entry data:", 
    sum(entry_combined$entry_source == "none", na.rm = TRUE), "\n")

cat("\n2. IML MERGE SUCCESS RATES\n")
cat("---------------------------\n")
cat("IML 1995 match rate:", 
    round(mean(!is.na(merged_cross_section$iml_95)) * 100, 1), "%\n")
cat("IML 2000 match rate:", 
    round(mean(!is.na(merged_cross_section$iml_00)) * 100, 1), "%\n")
cat("IML 2005 match rate:", 
    round(mean(!is.na(merged_cross_section$iml_05)) * 100, 1), "%\n")

cat("\n3. REGIONAL CLASSIFICATION\n")
cat("---------------------------\n")
cat("Localities with region assigned:", 
    sum(!is.na(merged_cross_section$region41)), "\n")
cat("Localities without region:", 
    sum(is.na(merged_cross_section$region41)), "\n")
cat("Number of unique regions:", 
    n_distinct(merged_cross_section$region41, na.rm = TRUE), "\n")

cat("\n4. FINAL DATASET SUMMARY\n")
cat("-------------------------\n")
cat("Cross-sectional dataset rows:", nrow(merged_final), "\n")
cat("Cross-sectional dataset columns:", ncol(merged_final), "\n")
cat("Panel dataset rows:", nrow(iml_panel_final), "\n")
cat("Panel dataset columns:", ncol(iml_panel_final), "\n")

cat("\n5. TREATMENT STATUS DISTRIBUTION\n")
cat("----------------------------------\n")
print(table(merged_final$program_phase))

cat("\n6. ENTRY YEAR DISTRIBUTION\n")
cat("----------------------------\n")
print(table(merged_final$entry_year))

sink()

message("\nStep 2: Dataset merging complete!")
message("Main dataset saved to: merged_entry_iml_rg.csv")
message("Panel dataset saved to: iml_panel_90_10.csv")
message("Next step: Run 03_create_rdd_datasets.R")