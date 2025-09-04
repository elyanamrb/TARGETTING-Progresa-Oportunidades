# =============================================================================
# Script: 03_assign_regions_and_geographic.R
# Purpose: Assign program regions to ENCASEH data and explore geographic patterns
#          - Map municipalities to program regions
#          - Handle missing region data, especially for Ventanilla
#          - Create geographic visualizations
#          - Document regional heterogeneity
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(haven)
library(scales)
library(viridis)
library(gridExtra)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/01-data-foundation")

# Create output directories
dir.create("outputs/regional_analysis", showWarnings = FALSE)

# Define consistent color scheme
color_scheme <- c(
  "Ventanilla" = "coral",
  "Rural Census" = "cornflowerblue",
  "Urban Census" = "steelblue"
)

# Theme for plots
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom"
  )

cat("========================================\n")
cat("REGIONAL ASSIGNMENT AND ANALYSIS\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD HARMONIZED DATA AND REGIONS FILE
# =============================================================================

cat("Loading harmonized data...\n")
encaseh_2000 <- readRDS("data/processed/encaseh_2000_harmonized.rds")
encaseh_2001 <- readRDS("data/processed/encaseh_2001_harmonized.rds")
encaseh_2002 <- readRDS("data/processed/encaseh_2002_harmonized.rds")
encaseh_2004 <- readRDS("data/processed/encaseh_2004_harmonized.rds")

# Load regions mapping file - adjust path as needed
cat("Loading regions mapping file...\n")
regions_map <- tryCatch({
  read_dta("data/raw/regiones41.dta")
}, error = function(e) {
  cat("WARNING: Could not load regiones41.dta\n")
  cat("Please ensure the file is in data/raw/\n")
  return(NULL)
})

# =============================================================================
# PART 2: EXPLORE REGIONS MAPPING
# =============================================================================

if (!is.null(regions_map)) {
  cat("\n========================================\n")
  cat("EXPLORING REGIONS MAPPING\n")
  cat("========================================\n")
  
  # Check structure
  cat("\nRegions mapping structure:\n")
  cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(regions_map), ncol(regions_map)))
  cat("Variables:", paste(names(regions_map), collapse = ", "), "\n")
  
  # Prepare regions mapping
  regions_map <- regions_map %>%
    mutate(
      cve_edo = as.numeric(as.character(cve_edo)),
      cve_mun = as.numeric(as.character(cve_mun)),
      region41 = as.numeric(as.character(region41)),
      muni_key = paste0(cve_edo, "_", cve_mun)
    )
  
  # Summary of regions
  n_regions <- length(unique(regions_map$region41))
  n_states <- length(unique(regions_map$cve_edo))
  n_municipalities <- nrow(regions_map)
  
  cat(sprintf("\nRegions mapping summary:\n"))
  cat(sprintf("  Number of regions: %d\n", n_regions))
  cat(sprintf("  Number of states: %d\n", n_states))
  cat(sprintf("  Number of municipalities: %d\n", n_municipalities))
  
  # Distribution of municipalities per region
  munis_per_region <- regions_map %>%
    count(region41) %>%
    summarise(
      mean = mean(n),
      median = median(n),
      min = min(n),
      max = max(n)
    )
  
  cat("\nMunicipalities per region:\n")
  cat(sprintf("  Mean: %.1f\n", munis_per_region$mean))
  cat(sprintf("  Median: %.0f\n", munis_per_region$median))
  cat(sprintf("  Range: %d - %d\n", munis_per_region$min, munis_per_region$max))
}

# =============================================================================
# PART 3: ASSIGN REGIONS TO ENCASEH DATA
# =============================================================================

cat("\n========================================\n")
cat("ASSIGNING REGIONS TO ENCASEH DATA\n")
cat("========================================\n")

assign_regions <- function(data, year, regions_df) {
  if (is.null(data) || is.null(regions_df)) {
    return(data)
  }
  
  cat(sprintf("\nAssigning regions to year %d...\n", year))
  
  # Check current region status
  n_with_region <- sum(!is.na(data$REGION))
  cat(sprintf("  Observations with existing REGION: %d (%.1f%%)\n", 
              n_with_region, n_with_region/nrow(data)*100))
  
  # Create municipality key
  data <- data %>%
    mutate(
      ENT_FED = as.numeric(as.character(ENT_FED)),
      MUNICI = as.numeric(as.character(MUNICI)),
      muni_key = paste0(ENT_FED, "_", MUNICI)
    )
  
  # Join with regions mapping
  data <- data %>%
    left_join(
      regions_df %>% select(muni_key, region41),
      by = "muni_key"
    )
  
  # Update REGION variable
  data$REGION_assigned <- data$region41
  
  # Use existing REGION if available, otherwise use assigned
  data$REGION_final <- ifelse(is.na(data$REGION), data$REGION_assigned, data$REGION)
  
  # Report assignment results
  n_assigned <- sum(!is.na(data$REGION_assigned) & is.na(data$REGION))
  n_final <- sum(!is.na(data$REGION_final))
  
  cat(sprintf("  Newly assigned regions: %d\n", n_assigned))
  cat(sprintf("  Total with regions: %d (%.1f%%)\n", 
              n_final, n_final/nrow(data)*100))
  
  # Check missing regions by targeting method
  if ("targeting_method" %in% names(data)) {
    missing_by_method <- data %>%
      group_by(targeting_method) %>%
      summarise(
        total = n(),
        missing_region = sum(is.na(REGION_final)),
        pct_missing = missing_region / total * 100,
        .groups = "drop"
      )
    
    cat("\n  Missing regions by targeting method:\n")
    print(missing_by_method)
  }
  
  # Clean up temporary columns
  data <- data %>%
    select(-muni_key, -region41, -REGION_assigned) %>%
    rename(REGION_original = REGION,
           REGION = REGION_final)
  
  return(data)
}

# Assign regions to each year
if (!is.null(regions_map)) {
  encaseh_2000 <- assign_regions(encaseh_2000, 2000, regions_map)
  encaseh_2001 <- assign_regions(encaseh_2001, 2001, regions_map)
  encaseh_2002 <- assign_regions(encaseh_2002, 2002, regions_map)
  encaseh_2004 <- assign_regions(encaseh_2004, 2004, regions_map)
}

# =============================================================================
# PART 4: ANALYZE REGIONAL PATTERNS
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING REGIONAL PATTERNS\n")
cat("========================================\n")

analyze_regional_patterns <- function(data, year) {
  if (is.null(data) || !"REGION" %in% names(data)) {
    return(NULL)
  }
  
  cat(sprintf("\nRegional analysis for year %d...\n", year))
  
  # Overall regional distribution
  regional_summary <- data %>%
    filter(!is.na(REGION)) %>%
    group_by(REGION) %>%
    summarise(
      n_households = n(),
      accept_rate = mean(accepted, na.rm = TRUE) * 100,
      mean_assets = mean(asset_count, na.rm = TRUE),
      mean_pun_inc = mean(pun_inc_clean, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = year)
  
  # If TIPO_CAP exists, analyze by targeting method
  if ("targeting_method" %in% names(data) && year >= 2001) {
    method_summary <- data %>%
      filter(!is.na(REGION) & !is.na(targeting_method)) %>%
      group_by(REGION, targeting_method) %>%
      summarise(
        n = n(),
        accept_rate = mean(accepted, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = targeting_method,
        values_from = c(n, accept_rate),
        values_fill = list(n = 0, accept_rate = NA)
      )
    
    # Calculate Ventanilla penetration for 2004
    if (year == 2004 && "n_Ventanilla" %in% names(method_summary)) {
      method_summary <- method_summary %>%
        mutate(
          total_n = rowSums(select(., starts_with("n_")), na.rm = TRUE),
          ventanilla_pct = n_Ventanilla / total_n * 100
        )
    }
    
    regional_summary <- left_join(regional_summary, method_summary, by = "REGION")
  }
  
  return(regional_summary)
}

# Analyze each year
regional_2000 <- analyze_regional_patterns(encaseh_2000, 2000)
regional_2001 <- analyze_regional_patterns(encaseh_2001, 2001)
regional_2002 <- analyze_regional_patterns(encaseh_2002, 2002)
regional_2004 <- analyze_regional_patterns(encaseh_2004, 2004)

# =============================================================================
# PART 5: CREATE REGIONAL VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("CREATING REGIONAL VISUALIZATIONS\n")
cat("========================================\n")

# 1. Ventanilla penetration by region (2004)
if (!is.null(regional_2004) && "ventanilla_pct" %in% names(regional_2004)) {
  
  # Top regions by Ventanilla usage
  top_ventanilla <- regional_2004 %>%
    filter(!is.na(ventanilla_pct) & ventanilla_pct > 0) %>%
    arrange(desc(ventanilla_pct)) %>%
    slice_head(n = 20)
  
  p1 <- ggplot(top_ventanilla, aes(x = reorder(factor(REGION), ventanilla_pct), 
                                   y = ventanilla_pct)) +
    geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
    coord_flip() +
    geom_text(aes(label = sprintf("%.1f%%", ventanilla_pct)), 
              hjust = -0.2, size = 3) +
    labs(
      title = "Top 20 Regions by Ventanilla Usage (2004)",
      subtitle = "Percentage of households using self-selection mechanism",
      x = "Region",
      y = "Ventanilla Usage (%)"
    ) +
    plot_theme +
    scale_y_continuous(limits = c(0, max(top_ventanilla$ventanilla_pct) * 1.1))
  
  ggsave("outputs/regional_analysis/ventanilla_penetration_2004.pdf", 
         p1, width = 10, height = 8)
  
  cat("✓ Created Ventanilla penetration plot\n")
}

# 2. Acceptance rates by region over time
if (!is.null(regional_2004)) {
  all_regional <- bind_rows(regional_2000, regional_2001, 
                            regional_2002, regional_2004)
  
  # Focus on regions with data in all years
  consistent_regions <- all_regional %>%
    group_by(REGION) %>%
    filter(n() == 4) %>%
    ungroup()
  
  if (nrow(consistent_regions) > 0) {
    # Select top 15 regions by size
    top_regions <- consistent_regions %>%
      group_by(REGION) %>%
      summarise(total_n = sum(n_households)) %>%
      arrange(desc(total_n)) %>%
      slice_head(n = 15) %>%
      pull(REGION)
    
    p2 <- consistent_regions %>%
      filter(REGION %in% top_regions) %>%
      ggplot(aes(x = factor(year), y = accept_rate, 
                 group = factor(REGION), color = factor(REGION))) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      facet_wrap(~REGION, scales = "free_y", ncol = 5) +
      labs(
        title = "Acceptance Rates by Region Over Time",
        subtitle = "Top 15 regions by sample size",
        x = "Year",
        y = "Acceptance Rate (%)"
      ) +
      plot_theme +
      theme(legend.position = "none")
    
    ggsave("outputs/regional_analysis/acceptance_rates_by_region_time.pdf", 
           p2, width = 12, height = 8)
    
    cat("✓ Created acceptance rates timeline\n")
  }
}

# 3. Regional heterogeneity in targeting methods (2004)
if (!is.null(encaseh_2004) && "targeting_method" %in% names(encaseh_2004)) {
  
  regional_composition <- encaseh_2004 %>%
    filter(!is.na(REGION) & !is.na(targeting_method)) %>%
    group_by(REGION) %>%
    mutate(region_total = n()) %>%
    group_by(REGION, targeting_method, region_total) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(percentage = n / region_total * 100) %>%
    filter(region_total >= 100)  # Focus on regions with sufficient data
  
  # Create stacked bar chart
  p3 <- ggplot(regional_composition, 
               aes(x = reorder(factor(REGION), -region_total), 
                   y = percentage, fill = targeting_method)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = color_scheme) +
    labs(
      title = "Targeting Method Composition by Region (2004)",
      subtitle = "Regions with at least 100 observations",
      x = "Region (ordered by sample size)",
      y = "Percentage (%)",
      fill = "Targeting Method"
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("outputs/regional_analysis/targeting_composition_by_region_2004.pdf", 
         p3, width = 12, height = 6)
  
  cat("✓ Created targeting method composition plot\n")
}

# 4. Heatmap of key indicators by region
if (!is.null(regional_2004)) {
  # Prepare data for heatmap
  heatmap_data <- regional_2004 %>%
    filter(!is.na(REGION)) %>%
    select(REGION, accept_rate, mean_assets, mean_pun_inc) %>%
    mutate(
      accept_rate_z = scale(accept_rate)[,1],
      assets_z = scale(mean_assets)[,1],
      pun_inc_z = scale(mean_pun_inc)[,1]
    ) %>%
    select(REGION, ends_with("_z")) %>%
    pivot_longer(cols = ends_with("_z"), 
                 names_to = "indicator", 
                 values_to = "z_score") %>%
    mutate(
      indicator = case_when(
        indicator == "accept_rate_z" ~ "Acceptance Rate",
        indicator == "assets_z" ~ "Asset Ownership",
        indicator == "pun_inc_z" ~ "Poverty Score",
        TRUE ~ indicator
      )
    )
  
  p4 <- ggplot(heatmap_data, aes(x = factor(REGION), y = indicator, fill = z_score)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred",
                         midpoint = 0, limits = c(-3, 3)) +
    labs(
      title = "Regional Characteristics Heatmap (2004)",
      subtitle = "Standardized values (z-scores)",
      x = "Region",
      y = "",
      fill = "Z-Score"
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 7))
  
  ggsave("outputs/regional_analysis/regional_characteristics_heatmap_2004.pdf", 
         p4, width = 14, height = 4)
  
  cat("✓ Created regional characteristics heatmap\n")
}

# =============================================================================
# PART 6: SAVE REGIONAL ANALYSIS RESULTS
# =============================================================================

cat("\n========================================\n")
cat("SAVING REGIONAL ANALYSIS RESULTS\n")
cat("========================================\n")

# Save updated datasets with regions
saveRDS(encaseh_2000, "data/processed/encaseh_2000_with_regions.rds")
saveRDS(encaseh_2001, "data/processed/encaseh_2001_with_regions.rds")
saveRDS(encaseh_2002, "data/processed/encaseh_2002_with_regions.rds")
saveRDS(encaseh_2004, "data/processed/encaseh_2004_with_regions.rds")

# Save regional summaries
all_regional_summary <- bind_rows(
  regional_2000,
  regional_2001,
  regional_2002,
  regional_2004
)

write.csv(all_regional_summary, 
          "outputs/regional_analysis/regional_summary_all_years.csv", 
          row.names = FALSE)

# Create key findings summary
if (!is.null(regional_2004) && "ventanilla_pct" %in% names(regional_2004)) {
  key_findings <- list(
    max_ventanilla_region = regional_2004 %>% 
      filter(ventanilla_pct == max(ventanilla_pct, na.rm = TRUE)) %>%
      pull(REGION),
    max_ventanilla_pct = max(regional_2004$ventanilla_pct, na.rm = TRUE),
    n_regions_with_ventanilla = sum(regional_2004$ventanilla_pct > 0, na.rm = TRUE),
    n_regions_no_ventanilla = sum(regional_2004$ventanilla_pct == 0, na.rm = TRUE)
  )
  
  cat("\n=== KEY REGIONAL FINDINGS ===\n")
  cat(sprintf("Maximum Ventanilla penetration: %.1f%% (Region %d)\n", 
              key_findings$max_ventanilla_pct, 
              key_findings$max_ventanilla_region))
  cat(sprintf("Regions with Ventanilla: %d\n", key_findings$n_regions_with_ventanilla))
  cat(sprintf("Regions without Ventanilla: %d\n", key_findings$n_regions_no_ventanilla))
}

cat("\n✓ Regional analysis complete!\n")
cat("✓ Updated datasets saved with regional assignments\n")
cat("✓ Visualizations saved to 'outputs/regional_analysis/'\n")
