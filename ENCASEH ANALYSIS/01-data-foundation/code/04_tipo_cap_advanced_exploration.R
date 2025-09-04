# =============================================================================
# Script: 04_tipo_cap_advanced_exploration.R
# Purpose: Advanced exploration and analysis of TIPO_CAP (targeting methods)
#          - Locality-level targeting mechanism analysis
#          - Regional heterogeneity in implementation
#          - Normalized RDD analysis across regions
#          - Temporal evolution patterns
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(haven)
library(scales)
library(gridExtra)
library(patchwork)
library(viridis)
library(rdrobust)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/01-data-foundation")

# Create output directories for advanced analysis
dir.create("outputs/tipo_cap_advanced", showWarnings = FALSE)
dir.create("outputs/tipo_cap_advanced/locality_level", showWarnings = FALSE)
dir.create("outputs/tipo_cap_advanced/regional_rdd", showWarnings = FALSE)
dir.create("outputs/tipo_cap_advanced/normalized_rdd", showWarnings = FALSE)

# Define consistent color scheme
color_scheme <- c(
  "Ventanilla Only" = "coral",
  "Rural Census Only" = "cornflowerblue",
  "Urban Census Only" = "#86c5da",
  "Mixed" = "darkgrey",
  "Ventanilla" = "coral",
  "Rural Census" = "cornflowerblue",
  "Urban Census" = "#86c5da"
)

# Plot theme
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "bottom"
  )

cat("========================================\n")
cat("ADVANCED TIPO_CAP EXPLORATION\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD DATA WITH REGIONS
# =============================================================================

cat("Loading data with regional assignments...\n")
encaseh_2000 <- readRDS("data/processed/encaseh_2000_with_regions.rds")
encaseh_2001 <- readRDS("data/processed/encaseh_2001_with_regions.rds")
encaseh_2002 <- readRDS("data/processed/encaseh_2002_with_regions.rds")
encaseh_2004 <- readRDS("data/processed/encaseh_2004_with_regions.rds")

# =============================================================================
# PART 2: LOCALITY-LEVEL TARGETING MECHANISM ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("LOCALITY-LEVEL TARGETING ANALYSIS\n")
cat("========================================\n")

analyze_targeting_by_locality <- function(data, year_val) {
  # Create a summary at locality level
  locality_summary <- data %>%
    filter(!is.na(TIPO_CAP)) %>%
    group_by(REGION, ENT_FED, MUNICI, LOCALI) %>%
    summarise(
      total_households = n(),
      n_ventanilla = sum(TIPO_CAP == 4, na.rm = TRUE),
      n_rural = sum(TIPO_CAP == 5, na.rm = TRUE),
      n_urban = sum(TIPO_CAP == 6, na.rm = TRUE),
      has_ventanilla = n_ventanilla > 0,
      has_rural = n_rural > 0,
      has_urban = n_urban > 0,
      n_targeting_types = sum(c(has_ventanilla, has_rural, has_urban)),
      targeting_category = case_when(
        n_targeting_types == 1 & has_ventanilla ~ "Ventanilla Only",
        n_targeting_types == 1 & has_rural ~ "Rural Census Only",
        n_targeting_types == 1 & has_urban ~ "Urban Census Only",
        n_targeting_types > 1 ~ "Mixed",
        TRUE ~ "Unknown"
      ),
      .groups = "drop"
    ) %>%
    mutate(year = year_val)
  
  # Create a region-level summary
  region_summary <- locality_summary %>%
    group_by(REGION, year) %>%
    summarise(
      total_localities = n(),
      ventanilla_only = sum(targeting_category == "Ventanilla Only"),
      rural_only = sum(targeting_category == "Rural Census Only"),
      urban_only = sum(targeting_category == "Urban Census Only"),
      mixed = sum(targeting_category == "Mixed"),
      pct_ventanilla_only = ventanilla_only / total_localities * 100,
      pct_rural_only = rural_only / total_localities * 100,
      pct_urban_only = urban_only / total_localities * 100,
      pct_mixed = mixed / total_localities * 100,
      .groups = "drop"
    ) %>%
    filter(total_localities > 0)
  
  return(list(
    locality_level = locality_summary,
    region_level = region_summary
  ))
}

# Apply to each year
results_2001 <- analyze_targeting_by_locality(encaseh_2001, 2001)
results_2002 <- analyze_targeting_by_locality(encaseh_2002, 2002)
results_2004 <- analyze_targeting_by_locality(encaseh_2004, 2004)

# Combine summaries
region_summary_all <- bind_rows(
  results_2001$region_level,
  results_2002$region_level,
  results_2004$region_level
)

locality_summary_all <- bind_rows(
  results_2001$locality_level,
  results_2002$locality_level,
  results_2004$locality_level
)

# Save summaries
write.csv(region_summary_all, 
          "outputs/tipo_cap_advanced/region_summary_all_years.csv", 
          row.names = FALSE)

write.csv(locality_summary_all, 
          "outputs/tipo_cap_advanced/locality_summary_all_years.csv", 
          row.names = FALSE)

cat("✓ Locality-level analysis complete\n")

# =============================================================================
# PART 3: REGIONAL HETEROGENEITY VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("CREATING REGIONAL HETEROGENEITY PLOTS\n")
cat("========================================\n")

# 1. Stacked bar charts by year
create_stacked_bar <- function(year_data, year_val) {
  ggplot(year_data, aes(x = factor(REGION), fill = targeting_category)) +
    geom_bar() +
    scale_fill_manual(values = color_scheme) +
    labs(
      title = paste("Targeting Mechanism Distribution by Region -", year_val),
      subtitle = "Number of localities by targeting type within each region",
      x = "Region",
      y = "Number of Localities",
      fill = "Targeting Type"
    ) +
    plot_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
}

for (year in c(2001, 2002, 2004)) {
  plot_data <- locality_summary_all %>% filter(year == !!year)
  p <- create_stacked_bar(plot_data, year)
  ggsave(sprintf("outputs/tipo_cap_advanced/locality_level/targeting_by_region_%d.pdf", year), 
         p, width = 12, height = 8)
}

cat("✓ Stacked bar charts created\n")

# 2. Regional heatmaps
prepare_heatmap_data <- function(region_data, year_val) {
  region_data %>%
    filter(year == year_val) %>%
    select(REGION, pct_ventanilla_only, pct_rural_only, pct_urban_only, pct_mixed) %>%
    pivot_longer(
      cols = starts_with("pct_"),
      names_to = "targeting_type",
      values_to = "percentage"
    ) %>%
    mutate(
      targeting_type = case_when(
        targeting_type == "pct_ventanilla_only" ~ "Ventanilla Only",
        targeting_type == "pct_rural_only" ~ "Rural Census Only",
        targeting_type == "pct_urban_only" ~ "Urban Census Only",
        targeting_type == "pct_mixed" ~ "Mixed"
      )
    )
}

for (year in c(2001, 2002, 2004)) {
  heatmap_data <- prepare_heatmap_data(region_summary_all, year)
  
  p <- ggplot(heatmap_data, aes(x = targeting_type, y = factor(REGION), fill = percentage)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", percentage)), color = "black", size = 2.5) +
    scale_fill_gradient(low = "white", high = "coral", limits = c(0, 100)) +
    labs(
      title = paste("Targeting Mechanism Distribution by Region -", year),
      subtitle = "Percentage of localities in each region by targeting type",
      x = "Targeting Type",
      y = "Region"
    ) +
    plot_theme +
    theme(axis.text.y = element_text(size = 7))
  
  ggsave(sprintf("outputs/tipo_cap_advanced/region_heatmap_%d.pdf", year), 
         p, width = 10, height = 14)
}

cat("✓ Regional heatmaps created\n")

# 3. Small multiples for evolution
small_multiples_data <- region_summary_all %>%
  pivot_longer(
    cols = c(pct_ventanilla_only, pct_rural_only, pct_urban_only, pct_mixed),
    names_to = "targeting_type",
    values_to = "percentage"
  ) %>%
  mutate(
    targeting_type = case_when(
      targeting_type == "pct_ventanilla_only" ~ "Ventanilla Only",
      targeting_type == "pct_rural_only" ~ "Rural Census Only",
      targeting_type == "pct_urban_only" ~ "Urban Census Only",
      targeting_type == "pct_mixed" ~ "Mixed"
    )
  )

p_small <- ggplot(small_multiples_data, 
                  aes(x = year, y = percentage, color = targeting_type, group = targeting_type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ REGION, scales = "free_y") +
  scale_color_manual(values = color_scheme) +
  labs(
    title = "Evolution of Targeting Mechanisms by Region (2001-2004)",
    subtitle = "Percentage of localities using each targeting method over time",
    x = "Year",
    y = "Percentage of Localities (%)",
    color = "Targeting Method"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  plot_theme +
  theme(
    strip.background = element_rect(fill = "lightgrey", color = NA),
    strip.text = element_text(face = "bold", size = 8)
  )

ggsave("outputs/tipo_cap_advanced/targeting_evolution_small_multiples.pdf", 
       p_small, width = 15, height = 12)

cat("✓ Small multiples plot created\n")

# =============================================================================
# PART 4: REGIONAL PUN_INC ANALYSIS BY TARGETING METHOD
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING PUN_INC BY REGION AND METHOD\n")
cat("========================================\n")

regional_pun_analysis <- function(data, year_val) {
  region_pun_summary <- data %>%
    filter(!is.na(TIPO_CAP), !is.na(PUN_INC)) %>%
    group_by(REGION, TIPO_CAP) %>%
    summarise(
      count = n(),
      mean_PUN_INC = mean(PUN_INC, na.rm = TRUE),
      median_PUN_INC = median(PUN_INC, na.rm = TRUE),
      sd_PUN_INC = sd(PUN_INC, na.rm = TRUE),
      q25 = quantile(PUN_INC, 0.25, na.rm = TRUE),
      q75 = quantile(PUN_INC, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      targeting_method = case_when(
        TIPO_CAP == 4 ~ "Ventanilla",
        TIPO_CAP == 5 ~ "Rural Census",
        TIPO_CAP == 6 ~ "Urban Census"
      ),
      year = year_val
    )
  
  return(region_pun_summary)
}

# Apply to each year
region_pun_2001 <- regional_pun_analysis(encaseh_2001, 2001)
region_pun_2002 <- regional_pun_analysis(encaseh_2002, 2002)
region_pun_2004 <- regional_pun_analysis(encaseh_2004, 2004)

region_pun_all <- bind_rows(region_pun_2001, region_pun_2002, region_pun_2004)

# Save PUN_INC analysis
write.csv(region_pun_all, 
          "outputs/tipo_cap_advanced/regional_pun_inc_by_method.csv", 
          row.names = FALSE)

# Create heatmap for 2004
region_pun_heatmap <- region_pun_2004 %>%
  filter(count >= 50) %>%
  ggplot(aes(x = targeting_method, y = factor(REGION), fill = mean_PUN_INC)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", mean_PUN_INC)), color = "black", size = 2.5) +
  scale_fill_gradient2(low = "cornflowerblue", mid = "white", high = "coral",
                       midpoint = median(region_pun_2004$mean_PUN_INC, na.rm = TRUE)) +
  labs(
    title = "Mean Poverty Scores by Region and Targeting Method (2004)",
    subtitle = "Higher values indicate greater poverty",
    x = "Targeting Method",
    y = "Region",
    fill = "Mean PUN_INC"
  ) +
  plot_theme +
  theme(axis.text.y = element_text(size = 7))

ggsave("outputs/tipo_cap_advanced/region_pun_inc_heatmap_2004.pdf", 
       region_pun_heatmap, width = 12, height = 14)

cat("✓ Regional PUN_INC analysis complete\n")

# =============================================================================
# PART 5: NORMALIZED RDD ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("RUNNING NORMALIZED RDD ANALYSIS\n")
cat("========================================\n")

# Function to find optimal cutoff for each region
find_regional_cutoffs <- function(data, year) {
  regions <- sort(unique(data$REGION[!is.na(data$REGION)]))
  cutoffs <- data.frame(REGION = regions, optimal_cutoff = NA)
  
  for (i in 1:nrow(cutoffs)) {
    region_num <- cutoffs$REGION[i]
    region_data <- data %>% filter(REGION == region_num)
    
    if (nrow(region_data) < 50) next
    
    region_data$outcome <- ifelse(region_data$accepted == 1, 1, 0)
    
    # Find optimal cutoff
    cutoff_range <- quantile(region_data$PUN_INC, 
                             probs = seq(0.2, 0.8, by = 0.01), 
                             na.rm = TRUE)
    max_effect <- 0
    optimal_cutoff <- NA
    
    for (cutoff in cutoff_range) {
      region_data$treatment <- ifelse(region_data$PUN_INC >= cutoff, 1, 0)
      tryCatch({
        model <- lm(outcome ~ treatment, data = region_data)
        coef <- abs(coef(model)["treatment"])
        if (!is.na(coef) && coef > max_effect) {
          max_effect <- coef
          optimal_cutoff <- cutoff
        }
      }, error = function(e) NULL)
    }
    
    cutoffs$optimal_cutoff[i] <- optimal_cutoff
  }
  
  cutoffs <- cutoffs %>% filter(!is.na(optimal_cutoff))
  return(cutoffs)
}

# Function to perform pooled normalized RDD
pooled_normalized_rdd <- function(data, cutoffs, year, bandwidth = NULL) {
  pooled_data <- data.frame()
  
  for (i in 1:nrow(cutoffs)) {
    region_num <- cutoffs$REGION[i]
    cutoff <- cutoffs$optimal_cutoff[i]
    
    region_data <- data %>% 
      filter(REGION == region_num) %>%
      mutate(
        PUN_INC_normalized = PUN_INC - cutoff,
        outcome = ifelse(accepted == 1, 1, 0),
        treatment = ifelse(PUN_INC >= cutoff, 1, 0)
      )
    
    pooled_data <- rbind(pooled_data, region_data)
  }
  
  if (!is.null(bandwidth)) {
    pooled_data <- pooled_data %>%
      filter(abs(PUN_INC_normalized) <= bandwidth)
  }
  
  # Fit pooled model
  pooled_model <- lm(outcome ~ treatment + PUN_INC_normalized + 
                       treatment:PUN_INC_normalized, data = pooled_data)
  model_summary <- summary(pooled_model)
  
  # Create plot
  pooled_plot <- ggplot(pooled_data, aes(x = PUN_INC_normalized, y = outcome)) +
    stat_summary_bin(bins = 30, fun = "mean", geom = "point", color = "blue", size = 2) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    geom_smooth(data = subset(pooled_data, PUN_INC_normalized < 0),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(pooled_data, PUN_INC_normalized >= 0),
                method = "lm", se = TRUE, color = "darkblue") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Pooled Normalized RDD Analysis",
      subtitle = paste("Year:", year, "- All Regions with Valid Cutoffs"),
      x = "Normalized PUN_INC Score (0 = Threshold)",
      y = "Program Acceptance Rate"
    ) +
    plot_theme
  
  return(list(
    year = year,
    n_obs = nrow(pooled_data),
    n_regions = nrow(cutoffs),
    rd_estimate = coef(pooled_model)["treatment"],
    std_error = model_summary$coefficients["treatment", "Std. Error"],
    p_value = model_summary$coefficients["treatment", "Pr(>|t|)"],
    plot = pooled_plot
  ))
}

# Run normalized RDD for each year
normalized_results <- list()

for (year in c(2001, 2002, 2004)) {
  cat(sprintf("\nProcessing year %d...\n", year))
  
  data <- switch(as.character(year),
                 "2001" = encaseh_2001,
                 "2002" = encaseh_2002,
                 "2004" = encaseh_2004)
  
  cutoffs <- find_regional_cutoffs(data, year)
  
  if (nrow(cutoffs) > 0) {
    result <- pooled_normalized_rdd(data, cutoffs, year, bandwidth = 1)
    normalized_results[[as.character(year)]] <- result
    
    # Save plot
    ggsave(sprintf("outputs/tipo_cap_advanced/normalized_rdd/pooled_rdd_%d.pdf", year),
           result$plot, width = 10, height = 6)
    
    cat(sprintf("  RD Estimate: %.3f (SE: %.3f, p: %.3f)\n", 
                result$rd_estimate, result$std_error, result$p_value))
  }
}

# Create comparison across years
if (length(normalized_results) > 0) {
  comparison_data <- map_df(names(normalized_results), function(yr) {
    res <- normalized_results[[yr]]
    tibble(
      year = as.numeric(yr),
      rd_estimate = res$rd_estimate,
      std_error = res$std_error,
      p_value = res$p_value,
      n_regions = res$n_regions
    )
  })
  
  p_comparison <- ggplot(comparison_data, aes(x = factor(year), y = rd_estimate, group = 1)) +
    geom_line(size = 1.2) +
    geom_point(size = 3, color = "blue") +
    geom_errorbar(aes(ymin = rd_estimate - std_error, 
                      ymax = rd_estimate + std_error), width = 0.2) +
    labs(
      title = "Change in Program Acceptance Threshold Effect Over Time",
      subtitle = "Normalized Pooled RDD Analysis",
      x = "Year",
      y = "RD Estimate (Treatment Effect at Threshold)"
    ) +
    plot_theme
  
  ggsave("outputs/tipo_cap_advanced/normalized_rdd/comparison_across_years.pdf",
         p_comparison, width = 10, height = 6)
  
  write.csv(comparison_data, 
            "outputs/tipo_cap_advanced/normalized_rdd/comparison_summary.csv",
            row.names = FALSE)
}

cat("✓ Normalized RDD analysis complete\n")

# =============================================================================
# PART 6: VENTANILLA PENETRATION ANALYSIS
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING VENTANILLA PENETRATION\n")
cat("========================================\n")

# Regional Ventanilla analysis for 2004
regional_ventanilla_2004 <- encaseh_2004 %>%
  group_by(REGION) %>%
  summarise(
    total_households = n(),
    ventanilla_count = sum(TIPO_CAP == 4, na.rm = TRUE),
    ventanilla_pct = ventanilla_count / total_households * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(ventanilla_pct))

# Top regions plot
top_ventanilla <- regional_ventanilla_2004 %>%
  filter(total_households >= 100) %>%
  slice_head(n = 20)

p_ventanilla <- ggplot(top_ventanilla, 
                       aes(x = reorder(factor(REGION), ventanilla_pct), 
                           y = ventanilla_pct)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", ventanilla_pct)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    title = "Ventanilla Implementation by Region (2004)",
    subtitle = "Top 20 regions - Percentage of households using self-selection",
    x = "Region",
    y = "Ventanilla Usage (%)"
  ) +
  scale_y_continuous(limits = c(0, max(top_ventanilla$ventanilla_pct) * 1.1)) +
  plot_theme

ggsave("outputs/tipo_cap_advanced/ventanilla_penetration_2004.pdf",
       p_ventanilla, width = 10, height = 8)

# Temporal heatmap
regional_ventanilla_all <- bind_rows(
  encaseh_2001 %>%
    group_by(REGION) %>%
    summarise(
      total = n(),
      ventanilla_pct = sum(TIPO_CAP == 4, na.rm = TRUE) / n() * 100,
      year = 2001,
      .groups = "drop"
    ),
  encaseh_2002 %>%
    group_by(REGION) %>%
    summarise(
      total = n(),
      ventanilla_pct = sum(TIPO_CAP == 4, na.rm = TRUE) / n() * 100,
      year = 2002,
      .groups = "drop"
    ),
  encaseh_2004 %>%
    group_by(REGION) %>%
    summarise(
      total = n(),
      ventanilla_pct = sum(TIPO_CAP == 4, na.rm = TRUE) / n() * 100,
      year = 2004,
      .groups = "drop"
    )
) %>%
  filter(total >= 100)

p_heatmap_temporal <- ggplot(regional_ventanilla_all,
                             aes(x = factor(year), 
                                 y = reorder(factor(REGION), ventanilla_pct), 
                                 fill = ventanilla_pct)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", ventanilla_pct)), 
            color = "black", size = 2) +
  scale_fill_gradient(low = "white", high = "coral") +
  labs(
    title = "Ventanilla Implementation Across Regions and Years",
    subtitle = "Percentage of households using self-selection",
    x = "Year",
    y = "Region",
    fill = "Ventanilla %"
  ) +
  plot_theme +
  theme(axis.text.y = element_text(size = 6))

ggsave("outputs/tipo_cap_advanced/ventanilla_temporal_heatmap.pdf",
       p_heatmap_temporal, width = 10, height = 12)

cat("✓ Ventanilla penetration analysis complete\n")

# =============================================================================
# SUMMARY REPORT
# =============================================================================

cat("\n========================================\n")
cat("ANALYSIS SUMMARY\n")
cat("========================================\n")

# Key findings
cat("\n=== KEY FINDINGS ===\n")

# 1. Ventanilla emergence
ventanilla_2004_total <- sum(encaseh_2004$TIPO_CAP == 4, na.rm = TRUE)
ventanilla_2004_pct <- ventanilla_2004_total / nrow(encaseh_2004) * 100

cat(sprintf("\n1. Ventanilla Implementation (2004):\n"))
cat(sprintf("   - Total households: %d (%.1f%%)\n", 
            ventanilla_2004_total, ventanilla_2004_pct))
cat(sprintf("   - Number of regions with Ventanilla: %d\n",
            sum(regional_ventanilla_2004$ventanilla_count > 0)))
cat(sprintf("   - Maximum penetration: %.1f%% (Region %d)\n",
            max(regional_ventanilla_2004$ventanilla_pct),
            regional_ventanilla_2004$REGION[1]))

# 2. Locality patterns
mixed_localities_2004 <- sum(results_2004$locality_level$targeting_category == "Mixed")
total_localities_2004 <- nrow(results_2004$locality_level)

cat(sprintf("\n2. Locality-Level Patterns (2004):\n"))
cat(sprintf("   - Total localities analyzed: %d\n", total_localities_2004))
cat(sprintf("   - Localities with mixed methods: %d (%.1f%%)\n",
            mixed_localities_2004, 
            mixed_localities_2004/total_localities_2004*100))

# 3. Regional heterogeneity
regions_high_ventanilla <- sum(regional_ventanilla_2004$ventanilla_pct > 20, na.rm = TRUE)

cat(sprintf("\n3. Regional Heterogeneity:\n"))
cat(sprintf("   - Regions with >20%% Ventanilla: %d\n", regions_high_ventanilla))
cat(sprintf("   - Standard deviation in Ventanilla %%: %.1f\n",
            sd(regional_ventanilla_2004$ventanilla_pct)))

cat("\n✓ Advanced TIPO_CAP exploration complete!\n")
cat("✓ All outputs saved to 'outputs/tipo_cap_advanced/'\n")
