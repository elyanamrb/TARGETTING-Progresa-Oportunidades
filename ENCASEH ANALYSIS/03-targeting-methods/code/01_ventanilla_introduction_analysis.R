# =============================================================================
# Script: 01_ventanilla_introduction_analysis.R
# Purpose: Analyze the introduction of Ventanilla (self-selection) in 2004
#          - Document emergence patterns
#          - Regional heterogeneity
#          - Comparison with census methods
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(gridExtra)
library(scales)
library(viridis)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/03-targeting-methods")

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/ventanilla_analysis", showWarnings = FALSE)
dir.create("outputs/ventanilla_analysis/regional", showWarnings = FALSE)

# Define color scheme
color_scheme <- c(
  "Ventanilla" = "coral",
  "Rural Census" = "cornflowerblue",
  "Urban Census" = "steelblue"
)

# Plot theme
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

cat("========================================\n")
cat("VENTANILLA INTRODUCTION ANALYSIS\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

cat("Loading data...\n")
encaseh_2001 <- readRDS("../01-data-foundation/data/processed/encaseh_2001_with_regions.rds")
encaseh_2002 <- readRDS("../01-data-foundation/data/processed/encaseh_2002_with_regions.rds")
encaseh_2004 <- readRDS("../01-data-foundation/data/processed/encaseh_2004_with_regions.rds")

# =============================================================================
# PART 2: DOCUMENT VENTANILLA EMERGENCE
# =============================================================================

cat("\n========================================\n")
cat("DOCUMENTING VENTANILLA EMERGENCE\n")
cat("========================================\n")

# Check Ventanilla presence by year
check_ventanilla <- function(data, year) {
  n_ventanilla <- sum(data$targeting_method == "Ventanilla", na.rm = TRUE)
  n_total <- nrow(data)
  pct_ventanilla <- n_ventanilla / n_total * 100
  
  cat(sprintf("Year %d: %d Ventanilla cases (%.2f%%)\n", 
              year, n_ventanilla, pct_ventanilla))
  
  return(data.frame(
    year = year,
    n_ventanilla = n_ventanilla,
    n_total = n_total,
    pct_ventanilla = pct_ventanilla
  ))
}

ventanilla_timeline <- bind_rows(
  check_ventanilla(encaseh_2001, 2001),
  check_ventanilla(encaseh_2002, 2002),
  check_ventanilla(encaseh_2004, 2004)
)

write.csv(ventanilla_timeline, 
          "outputs/ventanilla_analysis/ventanilla_timeline.csv", 
          row.names = FALSE)

# =============================================================================
# PART 3: ANALYZE 2004 VENTANILLA CHARACTERISTICS
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING VENTANILLA CHARACTERISTICS (2004)\n")
cat("========================================\n")

# Focus on 2004 when Ventanilla appears
data_2004 <- encaseh_2004 %>%
  filter(!is.na(targeting_method))

# Compare characteristics across targeting methods
method_comparison <- data_2004 %>%
  group_by(targeting_method) %>%
  summarise(
    n = n(),
    percentage = n() / nrow(data_2004) * 100,
    
    # Acceptance rates
    acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
    
    # Poverty scores
    mean_pun_inc = mean(pun_inc_clean, na.rm = TRUE),
    median_pun_inc = median(pun_inc_clean, na.rm = TRUE),
    sd_pun_inc = sd(pun_inc_clean, na.rm = TRUE),
    
    # Asset ownership
    mean_assets = mean(asset_count, na.rm = TRUE),
    
    # Regional presence
    n_regions = n_distinct(REGION, na.rm = TRUE),
    
    .groups = "drop"
  )

cat("\n=== Method Comparison ===\n")
print(method_comparison)

write.csv(method_comparison, 
          "outputs/ventanilla_analysis/method_comparison_2004.csv", 
          row.names = FALSE)

# =============================================================================
# PART 4: REGIONAL HETEROGENEITY
# =============================================================================

cat("\n========================================\n")
cat("ANALYZING REGIONAL HETEROGENEITY\n")
cat("========================================\n")

# Regional Ventanilla penetration
regional_ventanilla <- data_2004 %>%
  filter(!is.na(REGION)) %>%
  group_by(REGION) %>%
  summarise(
    total = n(),
    n_ventanilla = sum(targeting_method == "Ventanilla"),
    n_rural = sum(targeting_method == "Rural Census"),
    n_urban = sum(targeting_method == "Urban Census"),
    pct_ventanilla = n_ventanilla / total * 100,
    pct_rural = n_rural / total * 100,
    pct_urban = n_urban / total * 100,
    
    # Acceptance rates by method
    accept_ventanilla = mean(accepted[targeting_method == "Ventanilla"], na.rm = TRUE) * 100,
    accept_rural = mean(accepted[targeting_method == "Rural Census"], na.rm = TRUE) * 100,
    accept_urban = mean(accepted[targeting_method == "Urban Census"], na.rm = TRUE) * 100,
    
    .groups = "drop"
  ) %>%
  filter(total >= 50)  # Ensure sufficient sample

# Identify high Ventanilla regions
high_ventanilla <- regional_ventanilla %>%
  filter(pct_ventanilla > 20) %>%
  arrange(desc(pct_ventanilla))

cat(sprintf("\nRegions with >20%% Ventanilla: %d\n", nrow(high_ventanilla)))
if (nrow(high_ventanilla) > 0) {
  cat("Top 5 regions:\n")
  print(head(high_ventanilla[, c("REGION", "pct_ventanilla", "n_ventanilla")], 5))
}

write.csv(regional_ventanilla, 
          "outputs/ventanilla_analysis/regional_ventanilla_2004.csv", 
          row.names = FALSE)

# =============================================================================
# PART 5: CREATE VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("CREATING VISUALIZATIONS\n")
cat("========================================\n")

# 1. Targeting method composition pie chart
p_composition <- data_2004 %>%
  count(targeting_method) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ggplot(aes(x = "", y = percentage, fill = targeting_method)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), 
            position = position_stack(vjust = 0.5), size = 5) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = color_scheme) +
  labs(
    title = "Targeting Method Distribution (2004)",
    subtitle = sprintf("Total households: %s", comma(nrow(data_2004))),
    fill = "Method"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    legend.position = "bottom"
  )

ggsave("outputs/ventanilla_analysis/method_composition_2004.pdf", 
       p_composition, width = 8, height = 8)

# 2. Regional Ventanilla heatmap
p_regional <- regional_ventanilla %>%
  filter(n_ventanilla > 0) %>%
  arrange(desc(pct_ventanilla)) %>%
  slice_head(n = 30) %>%
  ggplot(aes(x = reorder(factor(REGION), pct_ventanilla), y = pct_ventanilla)) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", pct_ventanilla)), 
            hjust = -0.1, size = 3) +
  coord_flip() +
  labs(
    title = "Ventanilla Penetration by Region (2004)",
    subtitle = "Top 30 regions with Ventanilla presence",
    x = "Region",
    y = "Percentage Using Ventanilla (%)"
  ) +
  scale_y_continuous(limits = c(0, max(regional_ventanilla$pct_ventanilla) * 1.1)) +
  plot_theme

ggsave("outputs/ventanilla_analysis/regional_ventanilla_penetration.pdf", 
       p_regional, width = 10, height = 10)

# 3. Acceptance rates comparison
p_acceptance <- method_comparison %>%
  ggplot(aes(x = targeting_method, y = acceptance_rate, fill = targeting_method)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.1f%%", acceptance_rate)), 
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = color_scheme) +
  labs(
    title = "Acceptance Rates by Targeting Method (2004)",
    subtitle = "Percentage of households accepted into program",
    x = "Targeting Method",
    y = "Acceptance Rate (%)"
  ) +
  scale_y_continuous(limits = c(0, max(method_comparison$acceptance_rate) * 1.2)) +
  plot_theme +
  theme(legend.position = "none")

ggsave("outputs/ventanilla_analysis/acceptance_rates_comparison.pdf", 
       p_acceptance, width = 10, height = 6)

# 4. Poverty score distributions
p_pun_dist <- data_2004 %>%
  filter(!is.na(pun_inc_clean)) %>%
  ggplot(aes(x = pun_inc_clean, fill = targeting_method)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = color_scheme) +
  labs(
    title = "Poverty Score Distribution by Targeting Method (2004)",
    subtitle = "Kernel density estimation of PUN_INC scores",
    x = "Poverty Score (PUN_INC)",
    y = "Density",
    fill = "Method"
  ) +
  plot_theme

ggsave("outputs/ventanilla_analysis/poverty_score_distributions.pdf", 
       p_pun_dist, width = 10, height = 6)

# 5. Regional composition stacked bar
regional_long <- regional_ventanilla %>%
  select(REGION, pct_ventanilla, pct_rural, pct_urban) %>%
  pivot_longer(cols = starts_with("pct_"), 
               names_to = "method", 
               values_to = "percentage") %>%
  mutate(method = case_when(
    method == "pct_ventanilla" ~ "Ventanilla",
    method == "pct_rural" ~ "Rural Census",
    method == "pct_urban" ~ "Urban Census"
  ))

p_stacked <- regional_long %>%
  filter(REGION %in% head(regional_ventanilla$REGION[order(-regional_ventanilla$total)], 20)) %>%
  ggplot(aes(x = factor(REGION), y = percentage, fill = method)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = color_scheme) +
  labs(
    title = "Targeting Method Composition by Region (2004)",
    subtitle = "Top 20 regions by sample size",
    x = "Region",
    y = "Percentage (%)",
    fill = "Method"
  ) +
  plot_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("outputs/ventanilla_analysis/regional_composition_stacked.pdf", 
       p_stacked, width = 12, height = 8)

# =============================================================================
# PART 6: STATISTICAL TESTS
# =============================================================================

cat("\n========================================\n")
cat("STATISTICAL COMPARISONS\n")
cat("========================================\n")

# Test differences in acceptance rates
if (length(unique(data_2004$targeting_method)) > 1) {
  # Chi-square test for independence
  contingency_table <- table(data_2004$targeting_method, data_2004$accepted)
  chi_test <- chisq.test(contingency_table)
  
  cat("\n=== Chi-square Test: Targeting Method vs Acceptance ===\n")
  cat(sprintf("Chi-square statistic: %.2f\n", chi_test$statistic))
  cat(sprintf("p-value: %.4f\n", chi_test$p.value))
  cat(sprintf("Interpretation: %s\n", 
              ifelse(chi_test$p.value < 0.05, 
                     "Significant difference in acceptance rates",
                     "No significant difference in acceptance rates")))
  
  # ANOVA for poverty scores
  anova_result <- aov(pun_inc_clean ~ targeting_method, data = data_2004)
  anova_summary <- summary(anova_result)
  
  cat("\n=== ANOVA: Poverty Scores by Targeting Method ===\n")
  print(anova_summary)
}

# =============================================================================
# SUMMARY REPORT
# =============================================================================

cat("\n========================================\n")
cat("VENTANILLA INTRODUCTION SUMMARY\n")
cat("========================================\n")

cat("\n=== KEY FINDINGS ===\n")
cat(sprintf("1. Ventanilla first appears in 2004\n"))
cat(sprintf("   - %d households (%.1f%%)\n", 
            ventanilla_timeline$n_ventanilla[3],
            ventanilla_timeline$pct_ventanilla[3]))

cat(sprintf("\n2. Acceptance rate differences:\n"))
cat(sprintf("   - Ventanilla: %.1f%%\n", 
            method_comparison$acceptance_rate[method_comparison$targeting_method == "Ventanilla"]))
cat(sprintf("   - Rural Census: %.1f%%\n", 
            method_comparison$acceptance_rate[method_comparison$targeting_method == "Rural Census"]))
cat(sprintf("   - Urban Census: %.1f%%\n", 
            method_comparison$acceptance_rate[method_comparison$targeting_method == "Urban Census"]))

cat(sprintf("\n3. Regional concentration:\n"))
cat(sprintf("   - Present in %d regions\n", 
            sum(regional_ventanilla$n_ventanilla > 0)))
cat(sprintf("   - Highest penetration: %.1f%%\n", 
            max(regional_ventanilla$pct_ventanilla)))

cat("\n✓ Ventanilla introduction analysis complete!\n")
cat("✓ Results saved to 'outputs/ventanilla_analysis/'\n")

# Save workspace
save.image("outputs/ventanilla_analysis/ventanilla_workspace.RData")
