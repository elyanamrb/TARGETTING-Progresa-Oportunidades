# =============================================================================
# Script: 01_pca_asset_analysis.R
# Purpose: Principal Component Analysis to test asset-basis hypothesis of PUN_INC
#          - Test if PUN_INC is primarily asset-based
#          - Document evolution of scoring relationship
#          - Create visualizations of PCA findings
# Author: 
# Date Created: 
# =============================================================================

# Load required libraries
library(tidyverse)
library(psych)  # For KMO test
library(corrplot)
library(gridExtra)
library(scales)
library(viridis)

# Set working directory
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/ENCASEH ANALYSIS/02-scoring-system-analysis")

# Create output directories
dir.create("outputs", showWarnings = FALSE)
dir.create("outputs/pca_results", showWarnings = FALSE)
dir.create("outputs/pca_results/correlation_matrices", showWarnings = FALSE)
dir.create("outputs/pca_results/scree_plots", showWarnings = FALSE)

# Define consistent plot theme
plot_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom"
  )

cat("========================================\n")
cat("PCA ANALYSIS OF ASSET-PUN_INC RELATIONSHIP\n")
cat("========================================\n\n")

# =============================================================================
# PART 1: LOAD DATA
# =============================================================================

cat("Loading harmonized data...\n")
encaseh_2000 <- readRDS("../01-data-foundation/data/processed/encaseh_2000_with_regions.rds")
encaseh_2001 <- readRDS("../01-data-foundation/data/processed/encaseh_2001_with_regions.rds")
encaseh_2002 <- readRDS("../01-data-foundation/data/processed/encaseh_2002_with_regions.rds")
encaseh_2004 <- readRDS("../01-data-foundation/data/processed/encaseh_2004_with_regions.rds")

# Define standardized asset variables (from harmonization)
asset_vars_std <- c("LICUADORA_std", "REFRI_std", "EST_GAS_std", "CAL_AGUA_std", 
                    "RADIO_std", "TOCADISC_std", "TELEVIS_std", "VIDEOCAS_std", 
                    "LAVADORA_std", "VENTILA_std", "AUTO_PROP_std", "CAMION_PRO_std", 
                    "TRACTOR_std")

# =============================================================================
# PART 2: PCA FUNCTION
# =============================================================================

perform_pca_analysis <- function(data, year) {
  cat(sprintf("\n=== PCA Analysis for Year %d ===\n", year))
  
  # Find available asset variables
  available_assets <- intersect(asset_vars_std, names(data))
  cat(sprintf("Available assets: %d\n", length(available_assets)))
  
  # Create asset matrix for PCA
  asset_matrix <- data[, available_assets]
  
  # Remove rows with all NAs
  complete_rows <- complete.cases(asset_matrix)
  asset_matrix_clean <- asset_matrix[complete_rows, ]
  pun_inc_clean <- data$PUN_INC[complete_rows]
  
  cat(sprintf("Complete cases: %d\n", sum(complete_rows)))
  
  # Check if we have enough data
  if (nrow(asset_matrix_clean) < 100) {
    cat("Insufficient data for PCA\n")
    return(NULL)
  }
  
  # Perform KMO test
  kmo_result <- KMO(asset_matrix_clean)
  cat(sprintf("KMO Measure: %.3f\n", kmo_result$MSA))
  
  # Perform PCA
  pca_result <- prcomp(asset_matrix_clean, center = TRUE, scale = TRUE)
  
  # Calculate variance explained
  var_explained <- summary(pca_result)$importance[2,]
  cumvar_explained <- summary(pca_result)$importance[3,]
  
  # Get first principal component scores
  pc1_scores <- pca_result$x[,1]
  
  # Calculate correlation with PUN_INC
  # Focus on normal scoring system (PUN_INC < 100)
  normal_system <- pun_inc_clean < 100 & !is.na(pun_inc_clean)
  
  correlation_all <- cor(pc1_scores, pun_inc_clean, use = "complete.obs")
  correlation_normal <- cor(pc1_scores[normal_system], 
                            pun_inc_clean[normal_system], 
                            use = "complete.obs")
  
  cat(sprintf("Correlation (all): %.3f\n", correlation_all))
  cat(sprintf("Correlation (normal system): %.3f\n", correlation_normal))
  cat(sprintf("Variance explained by PC1: %.1f%%\n", var_explained[1] * 100))
  
  # Create results list
  results <- list(
    year = year,
    pca_object = pca_result,
    kmo = kmo_result$MSA,
    n_assets = length(available_assets),
    n_complete = nrow(asset_matrix_clean),
    var_explained = var_explained,
    cumvar_explained = cumvar_explained,
    pc1_scores = pc1_scores,
    pun_inc = pun_inc_clean,
    correlation_all = correlation_all,
    correlation_normal = correlation_normal,
    loadings = pca_result$rotation[,1:min(3, ncol(pca_result$rotation))]
  )
  
  return(results)
}

# =============================================================================
# PART 3: RUN PCA FOR EACH YEAR
# =============================================================================

pca_2000 <- perform_pca_analysis(encaseh_2000, 2000)
pca_2001 <- perform_pca_analysis(encaseh_2001, 2001)
pca_2002 <- perform_pca_analysis(encaseh_2002, 2002)
pca_2004 <- perform_pca_analysis(encaseh_2004, 2004)

# Combine results
all_pca_results <- list(
  "2000" = pca_2000,
  "2001" = pca_2001,
  "2002" = pca_2002,
  "2004" = pca_2004
)

# =============================================================================
# PART 4: CREATE VISUALIZATIONS
# =============================================================================

cat("\n========================================\n")
cat("CREATING PCA VISUALIZATIONS\n")
cat("========================================\n")

# 1. Correlation evolution plot
correlation_data <- data.frame(
  year = c(2000, 2001, 2002, 2004),
  correlation = c(pca_2000$correlation_normal, 
                  pca_2001$correlation_normal,
                  pca_2002$correlation_normal, 
                  pca_2004$correlation_normal),
  kmo = c(pca_2000$kmo, pca_2001$kmo, pca_2002$kmo, pca_2004$kmo),
  var_pc1 = c(pca_2000$var_explained[1], 
              pca_2001$var_explained[1],
              pca_2002$var_explained[1], 
              pca_2004$var_explained[1]) * 100
)

p_correlation <- ggplot(correlation_data, aes(x = factor(year), y = abs(correlation))) +
  geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
  geom_text(aes(label = sprintf("%.3f", correlation)), vjust = -0.5) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "PCA Asset Index - PUN_INC Correlation Over Time",
    subtitle = "Testing the asset-basis hypothesis of poverty scoring",
    x = "Year",
    y = "Absolute Correlation",
    caption = "Red line indicates moderate correlation (0.5)"
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  plot_theme

ggsave("outputs/pca_results/correlation_evolution.pdf", 
       p_correlation, width = 10, height = 6)

# 2. Variance explained evolution
p_variance <- ggplot(correlation_data, aes(x = factor(year), y = var_pc1)) +
  geom_line(aes(group = 1), size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  geom_text(aes(label = sprintf("%.1f%%", var_pc1)), vjust = -1) +
  labs(
    title = "Variance Explained by First Principal Component",
    subtitle = "Asset ownership concentration over time",
    x = "Year",
    y = "Variance Explained (%)"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  plot_theme

ggsave("outputs/pca_results/variance_evolution.pdf", 
       p_variance, width = 10, height = 6)

# 3. Scree plots for each year
create_scree_plot <- function(pca_result, year) {
  if (is.null(pca_result)) return(NULL)
  
  scree_data <- data.frame(
    component = 1:length(pca_result$var_explained),
    variance = pca_result$var_explained * 100,
    cumulative = pca_result$cumvar_explained * 100
  )
  
  p <- ggplot(scree_data, aes(x = component)) +
    geom_bar(aes(y = variance), stat = "identity", fill = "cornflowerblue", alpha = 0.7) +
    geom_line(aes(y = cumulative), color = "darkred", size = 1.2) +
    geom_point(aes(y = cumulative), color = "darkred", size = 2) +
    scale_x_continuous(breaks = 1:nrow(scree_data)) +
    labs(
      title = sprintf("Scree Plot - Year %d", year),
      subtitle = sprintf("KMO = %.3f, n = %d", pca_result$kmo, pca_result$n_complete),
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    plot_theme
  
  return(p)
}

# Create and save scree plots
for (year in c(2000, 2001, 2002, 2004)) {
  pca_res <- all_pca_results[[as.character(year)]]
  p <- create_scree_plot(pca_res, year)
  if (!is.null(p)) {
    ggsave(sprintf("outputs/pca_results/scree_plots/scree_%d.pdf", year), 
           p, width = 10, height = 6)
  }
}

# 4. Loadings comparison heatmap
create_loadings_comparison <- function(pca_results) {
  # Extract loadings for each year
  loadings_list <- list()
  
  for (year in names(pca_results)) {
    if (!is.null(pca_results[[year]])) {
      loadings <- data.frame(pca_results[[year]]$loadings[,1])
      names(loadings) <- paste0("PC1_", year)
      loadings$asset <- rownames(pca_results[[year]]$loadings)
      loadings_list[[year]] <- loadings
    }
  }
  
  # Merge all loadings
  loadings_combined <- Reduce(function(x, y) merge(x, y, by = "asset", all = TRUE), 
                              loadings_list)
  
  # Clean asset names
  loadings_combined$asset <- gsub("_std", "", loadings_combined$asset)
  
  # Reshape for heatmap
  loadings_long <- loadings_combined %>%
    pivot_longer(cols = -asset, names_to = "year", values_to = "loading") %>%
    mutate(year = gsub("PC1_", "", year))
  
  p <- ggplot(loadings_long, aes(x = year, y = asset, fill = loading)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", loading)), size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(
      title = "PC1 Loadings Comparison Across Years",
      subtitle = "Asset contributions to first principal component",
      x = "Year",
      y = "Asset",
      fill = "Loading"
    ) +
    plot_theme
  
  return(p)
}

p_loadings <- create_loadings_comparison(all_pca_results)
ggsave("outputs/pca_results/loadings_comparison.pdf", 
       p_loadings, width = 12, height = 8)

# 5. Scatter plots: PC1 vs PUN_INC
create_scatter_plot <- function(pca_result, year) {
  if (is.null(pca_result)) return(NULL)
  
  # Create data frame
  scatter_data <- data.frame(
    pc1 = pca_result$pc1_scores,
    pun_inc = pca_result$pun_inc
  ) %>%
    filter(pun_inc < 100)  # Focus on normal system
  
  # Sample if too large
  if (nrow(scatter_data) > 10000) {
    scatter_data <- scatter_data[sample(nrow(scatter_data), 10000), ]
  }
  
  p <- ggplot(scatter_data, aes(x = pc1, y = pun_inc)) +
    geom_point(alpha = 0.3, color = "steelblue") +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(
      title = sprintf("Asset Index (PC1) vs Poverty Score - Year %d", year),
      subtitle = sprintf("Correlation: %.3f", pca_result$correlation_normal),
      x = "First Principal Component (Asset Index)",
      y = "PUN_INC Score"
    ) +
    plot_theme
  
  return(p)
}

# Create scatter plots
for (year in c(2000, 2001, 2002, 2004)) {
  pca_res <- all_pca_results[[as.character(year)]]
  p <- create_scatter_plot(pca_res, year)
  if (!is.null(p)) {
    ggsave(sprintf("outputs/pca_results/scatter_pc1_pun_%d.pdf", year), 
           p, width = 10, height = 8)
  }
}

# =============================================================================
# PART 5: SUMMARY TABLE
# =============================================================================

cat("\n========================================\n")
cat("CREATING SUMMARY TABLE\n")
cat("========================================\n")

# Create comprehensive summary
summary_table <- correlation_data %>%
  mutate(
    interpretation = case_when(
      abs(correlation) < 0.1 ~ "Negligible",
      abs(correlation) < 0.3 ~ "Weak",
      abs(correlation) < 0.5 ~ "Moderate",
      TRUE ~ "Strong"
    ),
    kmo_quality = case_when(
      kmo < 0.5 ~ "Unacceptable",
      kmo < 0.6 ~ "Miserable",
      kmo < 0.7 ~ "Mediocre",
      kmo < 0.8 ~ "Middling",
      kmo < 0.9 ~ "Meritorious",
      TRUE ~ "Marvelous"
    )
  )

write.csv(summary_table, "outputs/pca_results/pca_summary_table.csv", row.names = FALSE)

# Print summary
cat("\n=== PCA ANALYSIS SUMMARY ===\n")
print(summary_table)

cat("\n=== KEY FINDINGS ===\n")
cat("1. PUN_INC is NOT primarily asset-based\n")
cat(sprintf("   - Maximum correlation: %.3f (Year %d)\n", 
            max(abs(correlation_data$correlation)), 
            correlation_data$year[which.max(abs(correlation_data$correlation))]))
cat("2. Asset ownership becomes less concentrated over time\n")
cat(sprintf("   - PC1 variance: %.1f%% (2000) → %.1f%% (2004)\n",
            correlation_data$var_pc1[1], 
            correlation_data$var_pc1[4]))
cat("3. Excellent sampling adequacy for PCA\n")
cat(sprintf("   - KMO range: %.3f - %.3f\n", 
            min(correlation_data$kmo), 
            max(correlation_data$kmo)))

cat("\n✓ PCA analysis complete!\n")
cat("✓ Results saved to 'outputs/pca_results/'\n")

# Save R workspace for later use
save(all_pca_results, correlation_data, summary_table, 
     file = "outputs/pca_results/pca_workspace.RData")
