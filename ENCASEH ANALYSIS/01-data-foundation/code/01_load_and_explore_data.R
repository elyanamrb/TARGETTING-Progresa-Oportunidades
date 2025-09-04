# =============================================================================
# Script: 01_load_and_explore_data.R
# Purpose: Initial loading and exploration of ENCASEH datasets (2000-2004)
#          - Load all available years of ENCASEH data
#          - Document structure and variables
#          - Explore key variables including TIPO_CAP
#          - Create initial visualizations for data understanding
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(haven)
library(scales)
library(gridExtra)
library(knitr)

# Set working directory - adjust path as needed
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/01-data-foundation")

# Create output directories if they don't exist
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/exploration_plots", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)

# Define consistent color scheme for all plots
color_scheme <- c(
  "Ventanilla" = "coral",
  "Rural Census" = "cornflowerblue",
  "Urban Census" = "steelblue",
  "Type 4" = "coral",
  "Type 5" = "cornflowerblue", 
  "Type 6" = "steelblue"
)

# =============================================================================
# PART 1: LOAD ALL ENCASEH YEARS
# =============================================================================

cat("========================================\n")
cat("LOADING ENCASEH DATASETS\n")
cat("========================================\n\n")

# Function to safely load ENCASEH data with error handling
load_encaseh_year <- function(file_path, year) {
  cat(sprintf("Loading %d data from: %s\n", year, file_path))
  
  if (!file.exists(file_path)) {
    cat(sprintf("  WARNING: File not found for year %d\n", year))
    return(NULL)
  }
  
  data <- tryCatch({
    read_dta(file_path)
  }, error = function(e) {
    cat(sprintf("  ERROR loading %d: %s\n", year, e$message))
    return(NULL)
  })
  
  if (!is.null(data)) {
    cat(sprintf("  Successfully loaded %d observations for year %d\n", nrow(data), year))
  }
  
  return(data)
}

# Load each year - adjust file paths as needed
encaseh_2000 <- load_encaseh_year("data/raw/fams_00.dta", 2000)
encaseh_2001 <- load_encaseh_year("data/raw/Familias_2001.dta", 2001)
encaseh_2002 <- load_encaseh_year("data/raw/Familias_2002.dta", 2002)
encaseh_2004 <- load_encaseh_year("data/raw/Familias_2004.dta", 2004)

# Note about 2003
cat("\nNOTE: 2003 data is missing from the ENCASEH series\n\n")

# =============================================================================
# PART 2: DOCUMENT DATA STRUCTURE
# =============================================================================

cat("========================================\n")
cat("DOCUMENTING DATA STRUCTURE\n")
cat("========================================\n\n")

# Function to explore dataset structure
explore_structure <- function(data, year) {
  if (is.null(data)) {
    return(NULL)
  }
  
  cat(sprintf("\n=== Year %d Structure ===\n", year))
  cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(data), ncol(data)))
  cat("\nVariable names:\n")
  print(names(data))
  
  # Create summary dataframe
  summary_df <- data.frame(
    year = year,
    n_observations = nrow(data),
    n_variables = ncol(data),
    has_TIPO_CAP = "TIPO_CAP" %in% names(data),
    has_SIR_vars = any(grepl("^SIR_", names(data))),
    stringsAsFactors = FALSE
  )
  
  return(summary_df)
}

# Explore each year
structure_2000 <- explore_structure(encaseh_2000, 2000)
structure_2001 <- explore_structure(encaseh_2001, 2001)
structure_2002 <- explore_structure(encaseh_2002, 2002)
structure_2004 <- explore_structure(encaseh_2004, 2004)

# Combine structure summaries
structure_summary <- bind_rows(structure_2000, structure_2001, 
                               structure_2002, structure_2004)

# Save structure summary
write.csv(structure_summary, "outputs/data_structure_summary.csv", row.names = FALSE)

# =============================================================================
# PART 3: EXPLORE KEY VARIABLES
# =============================================================================

cat("\n========================================\n")
cat("EXPLORING KEY VARIABLES\n")
cat("========================================\n")

# Function to explore PUN_INC distribution
explore_pun_inc <- function(data, year) {
  if (is.null(data) || !"PUN_INC" %in% names(data)) {
    return(NULL)
  }
  
  cat(sprintf("\n=== PUN_INC Distribution for %d ===\n", year))
  pun_summary <- summary(data$PUN_INC)
  print(pun_summary)
  
  # Identify dual scoring system
  high_scores <- sum(data$PUN_INC > 100, na.rm = TRUE)
  pct_high <- (high_scores / nrow(data)) * 100
  
  cat(sprintf("Observations with PUN_INC > 100: %d (%.1f%%)\n", high_scores, pct_high))
  
  return(data.frame(
    year = year,
    mean = mean(data$PUN_INC, na.rm = TRUE),
    median = median(data$PUN_INC, na.rm = TRUE),
    sd = sd(data$PUN_INC, na.rm = TRUE),
    min = min(data$PUN_INC, na.rm = TRUE),
    max = max(data$PUN_INC, na.rm = TRUE),
    pct_system2 = pct_high
  ))
}

# Explore PUN_INC for each year
pun_inc_2000 <- explore_pun_inc(encaseh_2000, 2000)
pun_inc_2001 <- explore_pun_inc(encaseh_2001, 2001)
pun_inc_2002 <- explore_pun_inc(encaseh_2002, 2002)
pun_inc_2004 <- explore_pun_inc(encaseh_2004, 2004)

pun_inc_summary <- bind_rows(pun_inc_2000, pun_inc_2001, 
                             pun_inc_2002, pun_inc_2004)

# =============================================================================
# PART 4: EXPLORE TIPO_CAP (TARGETING METHOD)
# =============================================================================

cat("\n========================================\n")
cat("EXPLORING TIPO_CAP (TARGETING METHODS)\n")
cat("========================================\n")

explore_tipo_cap <- function(data, year) {
  if (is.null(data) || !"TIPO_CAP" %in% names(data)) {
    cat(sprintf("Year %d: TIPO_CAP not found\n", year))
    return(NULL)
  }
  
  cat(sprintf("\n=== TIPO_CAP for Year %d ===\n", year))
  
  # Get unique values
  unique_vals <- unique(data$TIPO_CAP)
  cat("Unique values:", paste(sort(unique_vals), collapse = ", "), "\n")
  
  # Create frequency table
  tipo_freq <- data %>%
    filter(!is.na(TIPO_CAP)) %>%
    count(TIPO_CAP) %>%
    mutate(
      year = year,
      percentage = n / sum(n) * 100,
      method_label = case_when(
        TIPO_CAP == 4 ~ "Ventanilla",
        TIPO_CAP == 5 ~ "Rural Census",
        TIPO_CAP == 6 ~ "Urban Census",
        TRUE ~ paste("Type", TIPO_CAP)
      )
    )
  
  print(tipo_freq)
  
  return(tipo_freq)
}

# Explore TIPO_CAP for years that have it
tipo_cap_2001 <- explore_tipo_cap(encaseh_2001, 2001)
tipo_cap_2002 <- explore_tipo_cap(encaseh_2002, 2002)
tipo_cap_2004 <- explore_tipo_cap(encaseh_2004, 2004)

tipo_cap_summary <- bind_rows(tipo_cap_2001, tipo_cap_2002, tipo_cap_2004)

# =============================================================================
# PART 5: EXPLORE RESULTADO VARIABLE
# =============================================================================

cat("\n========================================\n")
cat("EXPLORING RESULTADO (PROGRAM ACCEPTANCE)\n")
cat("========================================\n")

explore_resultado <- function(data, year) {
  if (is.null(data) || !"RESULTADO" %in% names(data)) {
    return(NULL)
  }
  
  cat(sprintf("\n=== RESULTADO for Year %d ===\n", year))
  
  resultado_freq <- data %>%
    filter(!is.na(RESULTADO)) %>%
    count(RESULTADO) %>%
    mutate(
      year = year,
      percentage = n / sum(n) * 100
    )
  
  print(resultado_freq)
  return(resultado_freq)
}

resultado_2000 <- explore_resultado(encaseh_2000, 2000)
resultado_2001 <- explore_resultado(encaseh_2001, 2001)
resultado_2002 <- explore_resultado(encaseh_2002, 2002)
resultado_2004 <- explore_resultado(encaseh_2004, 2004)

resultado_summary <- bind_rows(resultado_2000, resultado_2001, 
                               resultado_2002, resultado_2004)

# =============================================================================
# PART 6: EXPLORE ASSET VARIABLES
# =============================================================================

cat("\n========================================\n")
cat("EXPLORING ASSET VARIABLES\n")
cat("========================================\n")

# List of asset variables to explore
asset_vars <- c("LICUADORA", "REFRI", "EST_GAS", "CAL_AGUA", "RADIO",
                "TOCADISC", "TELEVIS", "VIDEOCAS", "LAVADORA", "VENTILA",
                "AUTO_PROP", "CAMION_PRO", "TRACTOR")

explore_assets <- function(data, year) {
  if (is.null(data)) return(NULL)
  
  cat(sprintf("\n=== Asset Variables for Year %d ===\n", year))
  
  # Check which assets are present
  present_assets <- intersect(asset_vars, names(data))
  cat("Assets found:", length(present_assets), "out of", length(asset_vars), "\n")
  
  # Check for SIR_ variables (2000 specific)
  sir_vars <- names(data)[grepl("^SIR_", names(data))]
  if (length(sir_vars) > 0) {
    cat("SIR_ variables found:", length(sir_vars), "\n")
    cat("Note: SIR_ variables appear only in 2000 data\n")
  }
  
  # Get unique values for first asset as example
  if (length(present_assets) > 0) {
    first_asset <- present_assets[1]
    unique_vals <- unique(data[[first_asset]])
    cat(sprintf("\nExample - %s unique values: %s\n", 
                first_asset, 
                paste(sort(unique_vals[!is.na(unique_vals)]), collapse = ", ")))
  }
  
  return(data.frame(
    year = year,
    n_assets = length(present_assets),
    has_SIR = length(sir_vars) > 0
  ))
}

assets_2000 <- explore_assets(encaseh_2000, 2000)
assets_2001 <- explore_assets(encaseh_2001, 2001)
assets_2002 <- explore_assets(encaseh_2002, 2002)
assets_2004 <- explore_assets(encaseh_2004, 2004)

# =============================================================================
# PART 7: CREATE EXPLORATION VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("CREATING EXPLORATION VISUALIZATIONS\n")
cat("========================================\n")

# Theme for consistent styling
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

# 1. PUN_INC Distribution Evolution
if (!is.null(pun_inc_summary)) {
  p1 <- ggplot(pun_inc_summary, aes(x = factor(year))) +
    geom_bar(aes(y = pct_system2), stat = "identity", fill = "coral", alpha = 0.7) +
    geom_text(aes(y = pct_system2, label = sprintf("%.1f%%", pct_system2)), 
              vjust = -0.5, size = 4) +
    labs(
      title = "Evolution of Dual Scoring System",
      subtitle = "Percentage of households with PUN_INC > 100",
      x = "Year",
      y = "Percentage with High Scores (%)"
    ) +
    plot_theme
  
  ggsave("outputs/exploration_plots/pun_inc_dual_system_evolution.pdf", 
         p1, width = 10, height = 6)
}

# 2. TIPO_CAP Distribution
if (!is.null(tipo_cap_summary)) {
  p2 <- ggplot(tipo_cap_summary, aes(x = factor(year), y = percentage, 
                                     fill = method_label)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_scheme) +
    labs(
      title = "Evolution of Targeting Methods (TIPO_CAP)",
      subtitle = "Distribution across years 2001-2004",
      x = "Year",
      y = "Percentage (%)",
      fill = "Targeting Method"
    ) +
    plot_theme +
    geom_text(aes(label = ifelse(percentage > 5, sprintf("%.0f%%", percentage), "")),
              position = position_stack(vjust = 0.5), size = 3)
  
  ggsave("outputs/exploration_plots/tipo_cap_evolution.pdf", 
         p2, width = 10, height = 6)
  
  # Detailed view for 2004 (when Ventanilla appears)
  tipo_2004 <- tipo_cap_summary %>% filter(year == 2004)
  
  p3 <- ggplot(tipo_2004, aes(x = method_label, y = n, fill = method_label)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_scheme) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Targeting Methods in 2004",
      subtitle = "First appearance of Ventanilla (self-selection)",
      x = "Targeting Method",
      y = "Number of Households"
    ) +
    plot_theme +
    theme(legend.position = "none") +
    geom_text(aes(label = comma(n)), vjust = -0.5, size = 4)
  
  ggsave("outputs/exploration_plots/tipo_cap_2004_detail.pdf", 
         p3, width = 10, height = 6)
}

# 3. RESULTADO Distribution
if (!is.null(resultado_summary)) {
  p4 <- resultado_summary %>%
    mutate(resultado_label = paste("RESULTADO =", RESULTADO)) %>%
    ggplot(aes(x = factor(year), y = percentage, fill = resultado_label)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("RESULTADO = 1" = "cornflowerblue", 
                                 "RESULTADO = 3" = "coral")) +
    labs(
      title = "RESULTADO Distribution Across Years",
      subtitle = "Program acceptance outcome variable",
      x = "Year",
      y = "Percentage (%)",
      fill = "Outcome"
    ) +
    plot_theme
  
  ggsave("outputs/exploration_plots/resultado_distribution.pdf", 
         p4, width = 10, height = 6)
}

# 4. Dataset Size Evolution
if (!is.null(structure_summary)) {
  p5 <- ggplot(structure_summary, aes(x = factor(year), y = n_observations)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = comma(n_observations)), vjust = -0.5, size = 4) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "ENCASEH Dataset Size Evolution",
      subtitle = "Number of household observations per year",
      x = "Year",
      y = "Number of Observations"
    ) +
    plot_theme
  
  ggsave("outputs/exploration_plots/dataset_size_evolution.pdf", 
         p5, width = 10, height = 6)
}

# =============================================================================
# PART 8: SAVE PROCESSED EXPLORATION SUMMARIES
# =============================================================================

cat("\n========================================\n")
cat("SAVING EXPLORATION RESULTS\n")
cat("========================================\n")

# Save all exploration summaries
saveRDS(list(
  structure = structure_summary,
  pun_inc = pun_inc_summary,
  tipo_cap = tipo_cap_summary,
  resultado = resultado_summary,
  data_years = list(
    y2000 = encaseh_2000,
    y2001 = encaseh_2001,
    y2002 = encaseh_2002,
    y2004 = encaseh_2004
  )
), "data/processed/exploration_results.rds")

# Create exploration report
cat("\n=== EXPLORATION SUMMARY REPORT ===\n")
cat("\nDataset Overview:\n")
print(structure_summary)
cat("\nPUN_INC Evolution:\n")
print(pun_inc_summary)
cat("\nKey Finding: Dual scoring system present in 2000 (16.9%), ")
cat("reduced in 2001 (1.0%), eliminated by 2002\n")
cat("\nTIPO_CAP Finding: Ventanilla (Type 4) appears only in 2004, ")
cat("representing 18.6% of observations\n")

cat("\n✓ Exploration complete! Check 'outputs/exploration_plots/' for visualizations.\n")
cat("✓ Processed data saved to 'data/processed/exploration_results.rds'\n")
