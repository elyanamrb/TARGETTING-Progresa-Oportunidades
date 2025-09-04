# ============================================================================
# 01_regional_rdd_analysis.R
# ============================================================================
# Purpose: Regional RDD analysis for all 41 regions with plots
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# This script performs RDD analysis for each region and saves all plots
# ============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(rdrobust)
library(kableExtra)
library(stargazer)
library(gridExtra)

# ============================================================================
# DATA LOADING
# ============================================================================

# Load the processed data
setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.3-regional-rdd-2002")
merged_entry_iml_rg <- read.csv("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

cat("========================================\n")
cat("REGIONAL RDD ANALYSIS\n")
cat("========================================\n\n")
cat("Data loaded:", nrow(merged_entry_iml_rg), "observations\n")
cat("Regions:", length(unique(merged_entry_iml_rg$region41)), "\n\n")

# ============================================================================
# STEP 1: FIND OPTIMAL CUTOFFS FOR EACH REGION
# ============================================================================

analyze_region_complete <- function(data, region_num, post_year = 2002) {
  "
  Complete RDD analysis for a single region
  Matches the original analyze_region_complete function
  "
  
  # Filter for this region
  region_data <- data %>% 
    filter(region41 == region_num) %>%
    filter(!is.na(iml_95) & !is.na(entry_year_30))
  
  # Check minimum observations
  if (nrow(region_data) < 30) return(NULL)
  
  # Create post-year variable
  post_var <- paste0("post_", post_year)
  region_data[[post_var]] <- ifelse(region_data$entry_year_30 >= post_year, 1, 0)
  
  # Find optimal cutoff
  cutoff_range <- seq(-1.5, -0.5, by = 0.05)
  max_effect <- 0
  optimal_cutoff <- NA
  optimal_stats <- NULL
  
  for (cutoff in cutoff_range) {
    region_data$treatment <- ifelse(region_data$iml_95 >= cutoff, 1, 0)
    tryCatch({
      model <- lm(as.formula(paste(post_var, "~ treatment")), data = region_data)
      coef_val <- abs(coef(model)["treatment"])
      if (!is.na(coef_val) && coef_val > max_effect) {
        max_effect <- coef_val
        optimal_cutoff <- cutoff
        optimal_stats <- summary(model)$coefficients["treatment", ]
      }
    }, error = function(e) NULL)
  }
  
  if (is.na(optimal_cutoff)) return(NULL)
  
  # Create final model with optimal cutoff
  region_data$treatment <- ifelse(region_data$iml_95 >= optimal_cutoff, 1, 0)
  final_model <- lm(as.formula(paste(post_var, "~ treatment")), data = region_data)
  model_stats <- summary(final_model)$coefficients["treatment", ]
  
  # Create RDD plot
  rdd_plot <- ggplot(region_data, aes(x = iml_95, y = !!sym(post_var))) +
    # Binned scatter points - BLUE as in original
    stat_summary_bin(bins = 30, fun = "mean", geom = "point", color = "blue", size = 2) +
    # Vertical line at cutoff - RED DASHED as in original
    geom_vline(xintercept = optimal_cutoff, color = "red", linetype = "dashed") +
    # Fitted lines for each side - DARKBLUE as in original
    geom_smooth(data = subset(region_data, iml_95 < optimal_cutoff),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(region_data, iml_95 >= optimal_cutoff),
                method = "lm", se = TRUE, color = "darkblue") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = paste("Entry Timing RD Plot for Region", region_num),
         subtitle = paste("Optimal Cutoff:", round(optimal_cutoff, 2)),
         x = "Marginality Index (1995)",
         y = paste("Post", post_year, "Entry (0/1)")) +
    theme_minimal()
  
  # Create distribution plot
  dist_plot <- ggplot(region_data, aes(x = iml_95, fill = factor(!!sym(post_var)))) +
    geom_density(alpha = 0.5) +
    # Use PINK and TURQUOISE as in original
    scale_fill_manual(values = c("0" = "pink", "1" = "turquoise"),
                      labels = c(paste("Pre", post_year), paste("Post", post_year)),
                      name = paste("Post", post_year)) +
    labs(title = paste("IML Distribution by Entry Period - Region", region_num),
         subtitle = paste("N =", nrow(region_data)),
         x = "Marginality Index (1995)",
         y = "Density") +
    theme_minimal()
  
  return(list(
    region = region_num,
    optimal_cutoff = round(optimal_cutoff, 2),
    n_obs = nrow(region_data),
    rd_estimate = round(model_stats["Estimate"], 3),
    std_error = round(model_stats["Std. Error"], 3),
    p_value = model_stats["Pr(>|t|)"],
    rdd_plot = rdd_plot,
    dist_plot = dist_plot
  ))
}

# ============================================================================
# STEP 2: RUN ANALYSIS FOR ALL REGIONS AND SAVE PLOTS
# ============================================================================

# Create output directories matching original structure
output_dir <- "regional_plots_iml95_2002"
dir.create(output_dir, showWarnings = FALSE)

# Get all unique regions
regions <- sort(unique(merged_entry_iml_rg$region41))

# Initialize results storage
all_results <- list()
results_df <- data.frame()

cat("Analyzing regions and saving plots...\n")
cat("Output directory:", output_dir, "\n\n")

# Process each region
for (region in regions) {
  cat(sprintf("Processing Region %02d...", region))
  
  # Analyze region
  result <- analyze_region_complete(merged_entry_iml_rg, region, post_year = 2002)
  
  if (!is.null(result)) {
    # Store results
    all_results[[length(all_results) + 1]] <- result
    
    # Add to results dataframe
    results_df <- rbind(results_df, data.frame(
      region = result$region,
      optimal_cutoff = result$optimal_cutoff,
      n_obs = result$n_obs,
      rd_estimate = result$rd_estimate,
      std_error = result$std_error,
      p_value = ifelse(result$p_value < 0.001, "< 0.001", sprintf("%.3f", result$p_value))
    ))
    
    # Save plots with same naming convention as original
    tryCatch({
      # RDD plot
      ggsave(
        filename = file.path(output_dir, sprintf("region_%02d_rdd.pdf", result$region)),
        plot = result$rdd_plot,
        width = 10,
        height = 6
      )
      
      # Distribution plot
      ggsave(
        filename = file.path(output_dir, sprintf("region_%02d_dist.pdf", result$region)),
        plot = result$dist_plot,
        width = 10,
        height = 6
      )
      
      cat(" ✓ Plots saved\n")
      
    }, error = function(e) {
      cat(sprintf(" ✗ Failed to save plots: %s\n", e$message))
    })
  } else {
    cat(" - Skipped (insufficient data)\n")
  }
}

# Save summary CSV
write.csv(results_df, file.path(output_dir, "results_summary.csv"), row.names = FALSE)
cat("\n✓ Results summary saved to", file.path(output_dir, "results_summary.csv"), "\n")

# ============================================================================
# STEP 3: CREATE NORMALIZED REGIONAL RDD ANALYSIS
# ============================================================================

cat("\n========================================\n")
cat("NORMALIZED REGIONAL RDD ANALYSIS\n")
cat("========================================\n\n")

# Extract cutoffs for normalization
regional_cutoffs <- results_df %>%
  select(region, optimal_cutoff) %>%
  rename(cutoff = optimal_cutoff)

# Normalize data by regional cutoffs
normalized_data <- merged_entry_iml_rg %>%
  mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
  inner_join(regional_cutoffs, by = c("region41" = "region")) %>%
  mutate(normalized_iml = iml_95 - cutoff) %>%
  filter(!is.na(normalized_iml),
         !is.na(entry_year_30),
         entry_year_30 >= 1995 & entry_year_30 <= 2010)

# Run RDD for each region with normalized data
normalized_results <- data.frame()

for (region in unique(normalized_data$region41)) {
  region_data <- normalized_data %>% filter(region41 == region)
  
  if (nrow(region_data) >= 30) {
    # Run RDD with normalized running variable
    rd_result <- tryCatch({
      rdrobust(
        y = region_data$post_2002,
        x = region_data$normalized_iml,
        c = 0,
        kernel = "triangular",
        bwselect = "mserd"
      )
    }, error = function(e) NULL)
    
    if (!is.null(rd_result)) {
      normalized_results <- rbind(normalized_results, data.frame(
        region = region,
        rd_estimate = round(rd_result$coef[1], 3),
        std_error = round(rd_result$se[1], 3),
        p_value = round(rd_result$pv[1], 3),
        n_obs = nrow(region_data),
        bandwidth = round(rd_result$bws[1], 3)
      ))
    }
  }
}

# Save normalized results
write.csv(normalized_results, file.path(output_dir, "normalized_rd_results.csv"), row.names = FALSE)

# ============================================================================
# STEP 4: CREATE SUMMARY VISUALIZATION
# ============================================================================

cat("Creating summary visualizations...\n")

# Create coefficient plot for all regions
if (nrow(results_df) > 0) {
  # Add significance indicator
  results_df$significant <- ifelse(
    results_df$p_value == "< 0.001" | as.numeric(results_df$p_value) < 0.05,
    "Yes", "No"
  )
  
  # Create coefficient plot
  coef_plot <- ggplot(results_df, aes(x = factor(region), y = rd_estimate)) +
    geom_point(aes(color = significant), size = 3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("Yes" = "blue", "No" = "gray70"),
                       name = "Significant (p < 0.05)") +
    labs(title = "Regional RD Estimates",
         subtitle = "Effect of IML on Post-2002 Entry by Region",
         x = "Region",
         y = "RD Estimate") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(output_dir, "all_regions_coefficients.pdf"), 
         coef_plot, width = 12, height = 6)
}

# Create summary statistics table
summary_stats <- results_df %>%
  summarize(
    n_regions = n(),
    mean_effect = mean(rd_estimate),
    sd_effect = sd(rd_estimate),
    median_effect = median(rd_estimate),
    n_significant = sum(p_value == "< 0.001" | as.numeric(p_value) < 0.05, na.rm = TRUE),
    prop_significant = n_significant / n_regions
  )

cat("\n========================================\n")
cat("SUMMARY STATISTICS\n")
cat("========================================\n")
cat("Number of regions analyzed:", summary_stats$n_regions, "\n")
cat("Mean RD estimate:", round(summary_stats$mean_effect, 3), "\n")
cat("SD of RD estimates:", round(summary_stats$sd_effect, 3), "\n")
cat("Median RD estimate:", round(summary_stats$median_effect, 3), "\n")
cat("Significant regions (p < 0.05):", summary_stats$n_significant, 
    "(", round(summary_stats$prop_significant * 100, 1), "%)\n")

# ============================================================================
# STEP 5: CREATE LATEX TABLE FOR REGIONAL RESULTS
# ============================================================================

# Format results for LaTeX
results_latex <- results_df %>%
  mutate(
    significance = ifelse(p_value == "< 0.001" | as.numeric(p_value) < 0.01, "***",
                          ifelse(as.numeric(p_value) < 0.05, "**",
                                 ifelse(as.numeric(p_value) < 0.1, "*", ""))),
    rd_estimate_str = paste0(rd_estimate, significance)
  ) %>%
  select(region, optimal_cutoff, n_obs, rd_estimate_str, std_error, p_value)

# Create LaTeX table
latex_table <- paste0(
  "\\begin{table}[!htbp] \\centering\n",
  "\\caption{Regional RDD Results (2002 Threshold)}\n",
  "\\label{tab:regional_rdd}\n",
  "\\begin{tabular}{lccccc}\n",
  "\\hline\\hline\n",
  "Region & Cutoff & N & RD Estimate & Std. Error & p-value \\\\\n",
  "\\hline\n"
)

for (i in 1:nrow(results_latex)) {
  latex_table <- paste0(latex_table,
                        results_latex$region[i], " & ",
                        results_latex$optimal_cutoff[i], " & ",
                        results_latex$n_obs[i], " & ",
                        results_latex$rd_estimate_str[i], " & ",
                        results_latex$std_error[i], " & ",
                        results_latex$p_value[i], " \\\\\n"
  )
}

latex_table <- paste0(latex_table,
                      "\\hline\n",
                      "\\multicolumn{6}{l}{\\textit{Notes: *** p < 0.01, ** p < 0.05, * p < 0.1}} \\\\\n",
                      "\\hline\\hline\n",
                      "\\end{tabular}\n",
                      "\\end{table}"
)

writeLines(latex_table, file.path(output_dir, "regional_rdd_results.tex"))
cat("\n✓ LaTeX table saved to", file.path(output_dir, "regional_rdd_results.tex"), "\n")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

cat("\n========================================\n")
cat("REGIONAL ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Outputs saved in:", output_dir, "\n")
cat("- Individual region plots: region_XX_rdd.pdf and region_XX_dist.pdf\n")
cat("- Results summary: results_summary.csv\n")
cat("- Normalized results: normalized_rd_results.csv\n")
cat("- Coefficient plot: all_regions_coefficients.pdf\n")
cat("- LaTeX table: regional_rdd_results.tex\n")

# Save workspace for future reference
save.image(file = file.path(output_dir, "regional_analysis_workspace.RData"))
cat("\n✓ Complete workspace saved\n")