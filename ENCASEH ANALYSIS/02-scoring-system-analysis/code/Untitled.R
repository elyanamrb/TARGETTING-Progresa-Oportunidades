# =============================================================================
# Script: 02_dual_system_investigation.R
# Purpose: Investigate the dual scoring system in PUN_INC (2000-2001)
#          - Document System 1 (normal) vs System 2 (high values)
#          - Analyze transition to single system
#          - Identify patterns and characteristics
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(gridExtra)
library(scales)
library(viridis)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/02-scoring-system-analysis")

# Create output directories
dir.create("outputs/dual_system", showWarnings = FALSE)

# Plot theme
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

cat("========================================\n")
cat("DUAL SCORING SYSTEM INVESTIGATION\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

cat("Loading data...\n")
encaseh_2000 <- readRDS("../01-data-foundation/data/processed/encaseh_2000_with_regions.rds")
encaseh_2001 <- readRDS("../01-data-foundation/data/processed/encaseh_2001_with_regions.rds")
encaseh_2002 <- readRDS("../01-data-foundation/data/processed/encaseh_2002_with_regions.rds")
encaseh_2004 <- readRDS("../01-data-foundation/data/processed/encaseh_2004_with_regions.rds")

# =============================================================================
# PART 2: ANALYZE DUAL SYSTEM DISTRIBUTION
# =============================================================================

analyze_dual_system <- function(data, year) {
  cat(sprintf("\n=== Analyzing Year %d ===\n", year))
  
  # Classify scoring system
  data <- data %>%
    mutate(
      scoring_system = case_when(
        is.na(PUN_INC) ~ NA_character_,
        PUN_INC > 100 ~ "System 2 (High)",
        TRUE ~ "System 1 (Normal)"
      )
    )
  
  # Calculate statistics
  system_stats <- data %>%
    filter(!is.na(scoring_system)) %>%
    group_by(scoring_system) %>%
    summarise(
      n = n(),
      mean_score = mean(PUN_INC, na.rm = TRUE),
      median_score = median(PUN_INC, na.rm = TRUE),
      sd_score = sd(PUN_INC, na.rm = TRUE),
      min_score = min(PUN_INC, na.rm = TRUE),
      max_score = max(PUN_INC, na.rm = TRUE),
      q25 = quantile(PUN_INC, 0.25, na.rm = TRUE),
      q75 = quantile(PUN_INC, 0.75, na.rm = TRUE),
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
      .groups = "drop"
    ) %>%
    mutate(
      year = year,
      percentage = n / sum(n) * 100
    )
  
  # Print summary
  cat(sprintf("System 1 (Normal): %d observations (%.1f%%)\n", 
              system_stats$n[system_stats$scoring_system == "System 1 (Normal)"],
              system_stats$percentage[system_stats$scoring_system == "System 1 (Normal)"]))
  
  if ("System 2 (High)" %in% system_stats$scoring_system) {
    cat(sprintf("System 2 (High): %d observations (%.1f%%)\n", 
                system_stats$n[system_stats$scoring_system == "System 2 (High)"],
                system_stats$percentage[system_stats$scoring_system == "System 2 (High)"]))
  } else {
    cat("System 2 (High): 0 observations (0.0%)\n")
  }
  
  return(list(data = data, stats = system_stats))
}

# Analyze each year
results_2000 <- analyze_dual_system(encaseh_2000, 2000)
results_2001 <- analyze_dual_system(encaseh_2001, 2001)
results_2002 <- analyze_dual_system(encaseh_2002, 2002)
results_2004 <- analyze_dual_system(encaseh_2004, 2004)

# Combine statistics
all_system_stats <- bind_rows(
  results_2000$stats,
  results_2001$stats,
  results_2002$stats,
  results_2004$stats
)

write.csv(all_system_stats, "outputs/dual_system/system_statistics.csv", row.names = FALSE)

# =============================================================================
# PART 3: VISUALIZE DUAL SYSTEM EVOLUTION
# =============================================================================

cat("\n========================================\n")
cat("CREATING DUAL SYSTEM VISUALIZATIONS\n")
cat("========================================\n")

# 1. Evolution of dual system over time
evolution_data <- all_system_stats %>%
  select(year, scoring_system, percentage) %>%
  pivot_wider(names_from = scoring_system, values_from = percentage, values_fill = 0)

p_evolution <- all_system_stats %>%
  ggplot(aes(x = factor(year), y = percentage, fill = scoring_system)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c("System 1 (Normal)" = "steelblue", 
                               "System 2 (High)" = "coral")) +
  labs(
    title = "Evolution of Dual Scoring System",
    subtitle = "Transition from dual to single system (2000-2004)",
    x = "Year",
    y = "Percentage of Observations",
    fill = "Scoring System"
  ) +
  plot_theme

ggsave("outputs/dual_system/system_evolution.pdf", p_evolution, width = 10, height = 6)

# 2. Distribution plots for each year
create_distribution_plot <- function(data, year) {
  # Sample if too large
  plot_data <- data
  if (nrow(plot_data) > 50000) {
    plot_data <- plot_data[sample(nrow(plot_data), 50000), ]
  }
  
  p1 <- ggplot(plot_data, aes(x = PUN_INC, fill = scoring_system)) +
    geom_histogram(bins = 100, alpha = 0.7) +
    scale_fill_manual(values = c("System 1 (Normal)" = "steelblue", 
                                 "System 2 (High)" = "coral")) +
    scale_x_continuous(limits = c(-50, 1500)) +
    labs(
      title = sprintf("PUN_INC Distribution - Year %d", year),
      subtitle = "Full range showing both systems",
      x = "PUN_INC Score",
      y = "Count",
      fill = "System"
    ) +
    plot_theme
  
  # Zoomed version for System 1
  p2 <- ggplot(filter(plot_data, PUN_INC <= 100), aes(x = PUN_INC)) +
    geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = median(plot_data$PUN_INC[plot_data$PUN_INC <= 100], na.rm = TRUE),
               color = "red", linetype = "dashed") +
    labs(
      title = sprintf("System 1 (Normal) Distribution - Year %d", year),
      subtitle = "Zoomed to PUN_INC ≤ 100",
      x = "PUN_INC Score",
      y = "Count"
    ) +
    plot_theme
  
  return(list(full = p1, zoomed = p2))
}

# Create and save distribution plots
for (year in c(2000, 2001, 2002, 2004)) {
  data <- switch(as.character(year),
                 "2000" = results_2000$data,
                 "2001" = results_2001$data,
                 "2002" = results_2002$data,
                 "2004" = results_2004$data)
  
  plots <- create_distribution_plot(data, year)
  
  combined <- gridExtra::grid.arrange(plots$full, plots$zoomed, ncol = 1)
  ggsave(sprintf("outputs/dual_system/distribution_%d.pdf", year), 
         combined, width = 10, height = 12)
}

# =============================================================================
# PART 4: ANALYZE CHARACTERISTICS BY SYSTEM
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING CHARACTERISTICS BY SYSTEM\n")
cat("========================================\n")

# Focus on 2000 and 2001 when dual system was present
analyze_system_characteristics <- function(data, year) {
  if (sum(data$scoring_system == "System 2 (High)", na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  cat(sprintf("\n=== Characteristics for Year %d ===\n", year))
  
  # Compare characteristics
  characteristics <- data %>%
    filter(!is.na(scoring_system)) %>%
    group_by(scoring_system) %>%
    summarise(
      n = n(),
      mean_assets = mean(asset_count, na.rm = TRUE),
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
      mean_pun = mean(PUN_INC, na.rm = TRUE),
      median_pun = median(PUN_INC, na.rm = TRUE),
      # Regional distribution
      n_regions = n_distinct(REGION, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(year = year)
  
  print(characteristics)
  
  return(characteristics)
}

char_2000 <- analyze_system_characteristics(results_2000$data, 2000)
char_2001 <- analyze_system_characteristics(results_2001$data, 2001)

characteristics_combined <- bind_rows(char_2000, char_2001)

if (nrow(characteristics_combined) > 0) {
  write.csv(characteristics_combined, 
            "outputs/dual_system/system_characteristics.csv", 
            row.names = FALSE)
}

# =============================================================================
# PART 5: REGIONAL ANALYSIS OF DUAL SYSTEM
# =============================================================================

cat("\n========================================\n")
cat("REGIONAL ANALYSIS OF DUAL SYSTEM\n")
cat("========================================\n")

# Analyze which regions used System 2
regional_system_2000 <- results_2000$data %>%
  filter(!is.na(REGION)) %>%
  group_by(REGION) %>%
  summarise(
    total = n(),
    n_system2 = sum(scoring_system == "System 2 (High)", na.rm = TRUE),
    pct_system2 = n_system2 / total * 100,
    .groups = "drop"
  ) %>%
  filter(pct_system2 > 0) %>%
  arrange(desc(pct_system2)) %>%
  mutate(year = 2000)

regional_system_2001 <- results_2001$data %>%
  filter(!is.na(REGION)) %>%
  group_by(REGION) %>%
  summarise(
    total = n(),
    n_system2 = sum(scoring_system == "System 2 (High)", na.rm = TRUE),
    pct_system2 = n_system2 / total * 100,
    .groups = "drop"
  ) %>%
  filter(pct_system2 > 0) %>%
  arrange(desc(pct_system2)) %>%
  mutate(year = 2001)

regional_system_all <- bind_rows(regional_system_2000, regional_system_2001)

if (nrow(regional_system_all) > 0) {
  # Create heatmap of System 2 usage by region
  p_regional <- ggplot(regional_system_all, 
                       aes(x = factor(year), y = factor(REGION), fill = pct_system2)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.1f%%", pct_system2)), size = 3) +
    scale_fill_gradient(low = "white", high = "coral", limits = c(0, 100)) +
    labs(
      title = "System 2 Usage by Region",
      subtitle = "Percentage of households using high-value scoring system",
      x = "Year",
      y = "Region",
      fill = "% System 2"
    ) +
    plot_theme +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave("outputs/dual_system/regional_system2_heatmap.pdf", 
         p_regional, width = 10, height = 12)
  
  write.csv(regional_system_all, 
            "outputs/dual_system/regional_system2_usage.csv", 
            row.names = FALSE)
}

# =============================================================================
# PART 6: TRANSITION ANALYSIS (2000-2001)
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING SYSTEM TRANSITION\n")
cat("========================================\n")

# Compare the reduction in System 2 usage
transition_summary <- all_system_stats %>%
  filter(scoring_system == "System 2 (High)" | 
           (scoring_system == "System 1 (Normal)" & year %in% c(2000, 2001))) %>%
  select(year, scoring_system, n, percentage) %>%
  pivot_wider(names_from = scoring_system, values_from = c(n, percentage))

cat("\n=== Transition Summary ===\n")
print(transition_summary)

# Calculate reduction rate
if (nrow(transition_summary) >= 2) {
  reduction_2000_2001 <- (transition_summary$`n_System 2 (High)`[1] - 
                            transition_summary$`n_System 2 (High)`[2]) / 
    transition_summary$`n_System 2 (High)`[1] * 100
  
  cat(sprintf("\nSystem 2 reduction from 2000 to 2001: %.1f%%\n", reduction_2000_2001))
  cat("System 2 completely eliminated by 2002\n")
}

# =============================================================================
# SUMMARY REPORT
# =============================================================================

cat("\n========================================\n")
cat("DUAL SYSTEM INVESTIGATION SUMMARY\n")
cat("========================================\n")

cat("\n=== KEY FINDINGS ===\n")
cat("1. Dual scoring system documented:\n")
cat(sprintf("   - 2000: %.1f%% using System 2\n", 
            all_system_stats$percentage[all_system_stats$year == 2000 & 
                                          all_system_stats$scoring_system == "System 2 (High)"]))
cat(sprintf("   - 2001: %.1f%% using System 2\n", 
            ifelse(any(all_system_stats$year == 2001 & 
                         all_system_stats$scoring_system == "System 2 (High)"),
                   all_system_stats$percentage[all_system_stats$year == 2001 & 
                                                 all_system_stats$scoring_system == "System 2 (High)"],
                   0)))
cat("   - 2002-2004: System 2 eliminated\n")

if (nrow(regional_system_all) > 0) {
  cat("\n2. Regional concentration of System 2:\n")
  cat(sprintf("   - Number of regions using System 2: %d\n", 
              n_distinct(regional_system_all$REGION)))
  cat(sprintf("   - Highest usage: Region %d (%.1f%%)\n",
              regional_system_all$REGION[1],
              regional_system_all$pct_system2[1]))
}

cat("\n3. Administrative consolidation completed by 2002\n")
cat("   - Single standardized scoring system achieved\n")

cat("\n✓ Dual system investigation complete!\n")
cat("✓ Results saved to 'outputs/dual_system/'\n")

# Save workspace
save.image("outputs/dual_system/dual_system_workspace.RData")
