# 01_data_preparation.R
# Clean script for loading and preparing ENCASEH data

# Load required libraries
library(tidyverse)
library(haven)
library(data.table)

# Set working directory
setwd("/Users/elyanaramos/Desktop/TA:RA/RA/V4")


# Function to load and standardize ENCASEH data
load_encaseh_data <- function() {
  # Load all ENCASEH years
  encaseh00 <- read_dta("ENCASEH/fams_00.dta")
  encaseh01 <- read_dta("ENCASEH/Familias_2001.dta")
  encaseh02 <- read_dta("ENCASEH/Familias_2002.dta")
  encaseh04 <- read_dta("ENCASEH/Familias_2004.dta")
  
  # Load regions mapping
  regions <- read_dta("Padron familias/regiones41.dta")
  
  return(list(
    y2000 = encaseh00,
    y2001 = encaseh01,
    y2002 = encaseh02,
    y2004 = encaseh04,
    regions = regions
  ))
}

# Function to assign regions based on municipality codes
assign_regions <- function(data, regions_df) {
  # Prepare regions mapping
  regions_df <- regions_df %>%
    mutate(
      cve_edo = as.numeric(as.character(cve_edo)),
      cve_mun = as.numeric(as.character(cve_mun)),
      region41 = as.numeric(as.character(region41)),
      muni_key = paste0(cve_edo, "_", cve_mun)
    )
  
  # Assign regions to data
  data <- data %>%
    mutate(
      ENT_FED = as.numeric(as.character(ENT_FED)),
      MUNICI = as.numeric(as.character(MUNICI)),
      muni_key = paste0(ENT_FED, "_", MUNICI)
    ) %>%
    left_join(
      regions_df %>% select(muni_key, region41),
      by = "muni_key"
    ) %>%
    mutate(REGION = region41) %>%
    select(-region41, -muni_key)
  
  return(data)
}

# Function to create analysis variables
prepare_analysis_vars <- function(data) {
  # First create basic variables
  data <- data %>%
    mutate(
      # Binary outcome for program acceptance
      accepted = ifelse(RESULTADO == 1, 1, 0),
      
      # Identify scoring system
      scoring_system = ifelse(PUN_INC > 100, "System 2 (High)", "System 1 (Normal)")
    )
  
  # Add targeting method only if TIPO_CAP exists
  if ("TIPO_CAP" %in% names(data)) {
    data <- data %>%
      mutate(
        targeting_method = case_when(
          TIPO_CAP == 4 ~ "Ventanilla",
          TIPO_CAP == 5 ~ "Rural Census",
          TIPO_CAP == 6 ~ "Urban Census",
          TRUE ~ paste("Type", TIPO_CAP)
        )
      )
  } else {
    # Add NA column for consistency
    data$targeting_method <- NA_character_
  }
  
  return(data)
}

# Main data preparation pipeline
prepare_all_data <- function() {
  # Load raw data
  data_list <- load_encaseh_data()
  
  # Assign regions and prepare variables for each year
  data_list$y2000 <- assign_regions(data_list$y2000, data_list$regions) %>% 
    prepare_analysis_vars()
  
  data_list$y2001 <- assign_regions(data_list$y2001, data_list$regions) %>% 
    prepare_analysis_vars()
  
  data_list$y2002 <- assign_regions(data_list$y2002, data_list$regions) %>% 
    prepare_analysis_vars()
  
  data_list$y2004 <- assign_regions(data_list$y2004, data_list$regions) %>% 
    prepare_analysis_vars()
  
  # Print summary statistics
  cat("Data preparation complete. Summary:\n")
  cat("==================================\n")
  
  years <- c(2000, 2001, 2002, 2004)
  data_names <- c("y2000", "y2001", "y2002", "y2004")
  
  for (i in 1:length(years)) {
    d <- data_list[[data_names[i]]]
    cat(sprintf("\nYear %d:\n", years[i]))
    cat(sprintf("  Total observations: %d\n", nrow(d)))
    cat(sprintf("  With regions assigned: %d (%.1f%%)\n", 
                sum(!is.na(d$REGION)), 
                sum(!is.na(d$REGION))/nrow(d)*100))
    
    if ("TIPO_CAP" %in% names(d)) {
      cat("  Targeting methods:\n")
      tipo_table <- table(d$targeting_method, useNA = "ifany")
      for (j in 1:length(tipo_table)) {
        cat(sprintf("    %s: %d\n", names(tipo_table)[j], tipo_table[j]))
      }
    }
  }
  
  return(data_list)
}

# Run the preparation
prepared_data <- prepare_all_data()

# Save prepared data for use in other scripts
saveRDS(prepared_data, "prepared_encaseh_data.rds")

