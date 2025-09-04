# ============================================================================
# 01_national_rdd_analysis_2002.R
# ============================================================================
# Purpose: Simplified national RDD analysis matching original methodology
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
# ============================================================================

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(rdrobust)
library(kableExtra)
library(stargazer)

# ============================================================================
# DATA LOADING
# ============================================================================

setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.2-national-rdd-2002")
merged_entry_iml_rg <- read.csv("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

# Basic data check
cat("Data loaded:", nrow(merged_entry_iml_rg), "observations\n")
cat("Regions:", length(unique(merged_entry_iml_rg$region41)), "\n\n")

# ============================================================================
# STEP 1: FIND OPTIMAL CUTOFFS FOR EACH REGION (MATCHING ORIGINAL)
# ============================================================================

# This function matches your analyze_region_complete function
find_regional_cutoffs <- function(data, post_year = 2002) {
  
  regions <- sort(unique(data$region41))
  regional_cutoffs <- data.frame()
  
  for (region_num in regions) {
    
    # Filter for this region (matching your original filters)
    region_data <- data %>% 
      filter(region41 == region_num) %>%
      filter(!is.na(iml_95) & !is.na(entry_year_30))
    
    # Skip small regions
    if (nrow(region_data) < 30) next
    
    # Create post-year variable
    region_data$post_year <- ifelse(region_data$entry_year_30 >= post_year, 1, 0)
    
    # Find optimal cutoff (matching your exact range and step)
    cutoff_range <- seq(-1.5, -0.5, by = 0.05)
    max_effect <- 0
    optimal_cutoff <- NA
    
    for (cutoff in cutoff_range) {
      region_data$treatment <- ifelse(region_data$iml_95 >= cutoff, 1, 0)
      
      tryCatch({
        model <- lm(post_year ~ treatment, data = region_data)
        coef_val <- abs(coef(model)["treatment"])
        
        if (!is.na(coef_val) && coef_val > max_effect) {
          max_effect <- coef_val
          optimal_cutoff <- cutoff
        }
      }, error = function(e) NULL)
    }
    
    if (!is.na(optimal_cutoff)) {
      regional_cutoffs <- rbind(regional_cutoffs,
                                data.frame(region = region_num,
                                           cutoff = optimal_cutoff))
    }
  }
  
  cat("Found cutoffs for", nrow(regional_cutoffs), "regions\n")
  return(regional_cutoffs)
}

# Find the regional cutoffs
cat("Finding optimal regional cutoffs...\n")
regional_cutoffs <- find_regional_cutoffs(merged_entry_iml_rg, post_year = 2002)

# ============================================================================
# STEP 2: NORMALIZE DATA BY REGIONAL CUTOFFS (MATCHING ORIGINAL)
# ============================================================================

# This matches your normalized_national_rdd function approach
cat("\nNormalizing by regional cutoffs...\n")

# Create post_2002 variable
merged_entry_iml_rg$post_2002 <- ifelse(merged_entry_iml_rg$entry_year_30 >= 2002, 1, 0)

# Join with regional cutoffs and normalize
normalized_data <- merged_entry_iml_rg %>%
  inner_join(regional_cutoffs, by = c("region41" = "region")) %>%
  mutate(normalized_iml = iml_95 - cutoff) %>%
  filter(!is.na(normalized_iml),
         !is.na(entry_year_30),
         entry_year_30 >= 1995 & entry_year_30 <= 2010)  # Matching your filter

cat("Normalized data:", nrow(normalized_data), "observations\n")
cat("Normalized IML range:", round(min(normalized_data$normalized_iml), 3), 
    "to", round(max(normalized_data$normalized_iml), 3), "\n\n")

# ============================================================================
# STEP 3: RUN NATIONAL RDD 
# ============================================================================

cat("Running national RDD with normalized data...\n\n")

# Run RDD with same specifications as your original
rd_result <- rdrobust(
  y = normalized_data$post_2002,
  x = normalized_data$normalized_iml,
  c = 0,  # Cutoff at 0 after normalization
  kernel = "triangular",
  bwselect = "mserd",
  all = TRUE
)

# Print summary matching your original output format
cat("=== SHARP RD ESTIMATES ===\n")
cat("Conventional estimate:", round(rd_result$coef[1], 3), "\n")
cat("Robust estimate:", round(rd_result$coef[2], 3), "\n")
cat("Std. Error (conv):", round(rd_result$se[1], 3), "\n")
cat("Std. Error (robust):", round(rd_result$se[2], 3), "\n")
cat("z-statistic (conv):", round(rd_result$z[1], 3), "\n")
cat("z-statistic (robust):", round(rd_result$z[2], 3), "\n")
cat("P-value (conv):", round(rd_result$pv[1], 4), "\n")
cat("P-value (robust):", round(rd_result$pv[2], 4), "\n")
cat("\nNumber of observations:", rd_result$N[1] + rd_result$N[2], "\n")
cat("Bandwidth:", round(rd_result$bws[1], 4), "\n")
cat("Observations left of cutoff:", rd_result$N[1], "\n")
cat("Observations right of cutoff:", rd_result$N[2], "\n")
cat("Effective N (left):", rd_result$N_h[1], "\n")
cat("Effective N (right):", rd_result$N_h[2], "\n\n")

# ============================================================================
# STEP 4: CREATE VISUALIZATIONS 
# ============================================================================

# Create binned scatter plot for RDD
create_rdd_plot <- function(data, rd_result) {
  
  # Create bins (matching your approach with stat_summary_bin)
  bin_data <- data %>%
    mutate(bin = cut(normalized_iml, breaks = 30, include.lowest = TRUE)) %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(normalized_iml),
      mean_y = mean(post_2002),
      n = n(),
      .groups = 'drop'
    ) %>%
    filter(n >= 5)
  
  # Create plot with your exact specifications
  p <- ggplot(bin_data, aes(x = mean_x, y = mean_y)) +
    # Blue points as specified
    geom_point(color = "blue", size = 2) +
    # Red dashed vertical line at cutoff
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    # Fitted lines for each side
    geom_smooth(data = subset(bin_data, mean_x < 0),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(bin_data, mean_x >= 0),
                method = "lm", se = TRUE, color = "darkblue") +
    # Scale and labels
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "National RDD: Program Entry Timing",
         subtitle = paste0("2002 Cutoff | Effect = ", round(rd_result$coef[1], 3)),
         x = "Normalized Marginality Index (1995)",
         y = "Post-2002 Entry (0/1)") +
    theme_minimal()
  
  return(p)
}

# Create distribution plot
create_distribution_plot <- function(data) {
  
  data$entry_period <- factor(
    ifelse(data$entry_year_30 >= 2002, "Post-2002", "Pre-2002"),
    levels = c("Pre-2002", "Post-2002")
  )
  
  p <- ggplot(data, aes(x = normalized_iml, fill = entry_period)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c("Pre-2002" = "coral", "Post-2002" = "turquoise"),
                      name = "Entry Period") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(title = "IML Distribution by Entry Period",
         subtitle = paste("N =", nrow(data)),
         x = "Normalized Marginality Index (1995)",
         y = "Density") +
    theme_minimal() +
    theme(legend.position = "top")
  
  return(p)
}

# Create and save plots
cat("Creating visualizations...\n")

# Create output directories if they don't exist
dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE)
dir.create("output/tables", showWarnings = FALSE)

# Generate plots
rdd_plot <- create_rdd_plot(normalized_data, rd_result)
dist_plot <- create_distribution_plot(normalized_data)

# Save plots
ggsave("output/figures/national_rdd_2002.pdf", rdd_plot, width = 10, height = 6)
ggsave("output/figures/national_rdd_2002.png", rdd_plot, width = 10, height = 6, dpi = 300)
ggsave("output/figures/national_distribution_2002.pdf", dist_plot, width = 10, height = 6)
ggsave("output/figures/national_distribution_2002.png", dist_plot, width = 10, height = 6, dpi = 300)

cat("Plots saved to output/figures/\n\n")

# ============================================================================
# STEP 5: CREATE LATEX TABLE
# ============================================================================

# Create results table
results_df <- data.frame(
  Method = c("Conventional", "Robust"),
  Coefficient = c(round(rd_result$coef[1], 3), 
                  ifelse(is.na(rd_result$coef[2]), "-", round(rd_result$coef[2], 3))),
  Std_Error = c(round(rd_result$se[1], 3), 
                ifelse(is.na(rd_result$se[2]), "-", round(rd_result$se[2], 3))),
  z_statistic = c(round(rd_result$z[1], 3), 
                  ifelse(is.na(rd_result$z[2]), "-", round(rd_result$z[2], 3))),
  p_value = c(ifelse(rd_result$pv[1] < 0.001, "< 0.001", round(rd_result$pv[1], 3)),
              ifelse(is.na(rd_result$pv[2]), "-", 
                     ifelse(rd_result$pv[2] < 0.001, "< 0.001", round(rd_result$pv[2], 3))))
)

# Save as CSV
write.csv(results_df, "output/tables/sharp_rd_estimates.csv", row.names = FALSE)

# Create LaTeX table manually
latex_table <- paste0(
  "\\begin{table}[!htbp] \\centering\n",
  "  \\caption{Table 1: Sharp RD Estimates Using Local Polynomial Regression}\n",
  "  \\label{tab:sharp_rd}\n",
  "\\begin{tabular}{lcccc}\n",
  "\\hline\\hline\n",
  "Method & Coefficient & Std. Error & z-statistic & p-value \\\\\n",
  "\\hline\n",
  "Conventional & ", results_df$Coefficient[1], " & ", results_df$Std_Error[1], 
  " & ", results_df$z_statistic[1], " & ", results_df$p_value[1], " \\\\\n",
  "Robust & ", results_df$Coefficient[2], " & ", results_df$Std_Error[2], 
  " & ", results_df$z_statistic[2], " & ", results_df$p_value[2], " \\\\\n",
  "\\hline\n",
  "\\multicolumn{5}{l}{\\textit{Notes:}} \\\\\n",
  "\\multicolumn{5}{l}{Number of Observations: ", format(rd_result$N[1] + rd_result$N[2], big.mark = ","), "} \\\\\n",
  "\\multicolumn{5}{l}{Bandwidth type: mserd} \\\\\n",
  "\\multicolumn{5}{l}{Kernel: Triangular} \\\\\n",
  "\\multicolumn{5}{l}{Number of Observations (left, right): ", 
  format(rd_result$N[1], big.mark = ","), ", ", format(rd_result$N[2], big.mark = ","), "} \\\\\n",
  "\\multicolumn{5}{l}{Effective Number of Observations (left, right): ", 
  format(rd_result$N_h[1], big.mark = ","), ", ", format(rd_result$N_h[2], big.mark = ","), "} \\\\\n",
  "\\hline\\hline\n",
  "\\end{tabular}\n",
  "\\end{table}"
)

writeLines(latex_table, "output/tables/sharp_rd_estimates.tex")
cat("LaTeX table saved to output/tables/sharp_rd_estimates.tex\n")

