# ============================================================================
# 01_yearly_rdd_evolution.R
# ============================================================================
# Purpose: Complete temporal analysis of RD estimates across threshold years
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# This script performs multi-year RDD analysis to understand how treatment
# effects evolve when using different entry year thresholds (1998-2006)
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rdrobust)
library(gridExtra)
library(grid)
library(kableExtra)
library(stargazer)
library(cowplot)

# ============================================================================
# SETUP
# ============================================================================

setwd("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.4-temporal-analysis")
merged_entry_iml_rg <- read.csv("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")


cat("========================================\n")
cat("TEMPORAL ANALYSIS - RD EVOLUTION\n")
cat("Issue 1.4: Multi-Year Analysis\n")
cat("========================================\n\n")

# Create output directories
dir.create("multi_year_analysis", showWarnings = FALSE)
dir.create("multi_year_analysis/regional_plots", showWarnings = FALSE)
dir.create("region_graphs", showWarnings = FALSE)

# ============================================================================
# CORE ANALYSIS FUNCTIONS 
# ============================================================================

# Function to find optimal cutoff for a specific region and year
analyze_region_rdd <- function(data, region_num = NULL, threshold_year, 
                               cutoff_range = seq(-1.5, -0.5, by = 0.05)) {
  # Filter by region if specified
  if (!is.null(region_num)) {
    data <- data %>% filter(region41 == region_num)
  }
  
  if (nrow(data) < 30) return(NULL)
  
  # Create outcome variable
  data <- data %>%
    mutate(post_year = ifelse(entry_year_30 >= threshold_year, 1, 0))
  
  # Find optimal cutoff
  max_effect <- -Inf
  optimal_cutoff <- NA
  optimal_result <- NULL
  
  for (cutoff in cutoff_range) {
    tryCatch({
      rd_result <- rdrobust(
        y = data$post_year,
        x = data$iml_95,
        c = cutoff,
        kernel = "triangular",
        bwselect = "mserd"
      )
      
      effect_size <- abs(rd_result$coef[1])
      
      if (!is.na(effect_size) && effect_size > max_effect) {
        max_effect <- effect_size
        optimal_cutoff <- cutoff
        optimal_result <- rd_result
      }
    }, error = function(e) NULL)
  }
  
  if (is.na(optimal_cutoff)) return(NULL)
  
  # Create plot data
  plot_data <- data %>%
    mutate(bin = cut(iml_95, breaks = 50))
  
  bin_means <- plot_data %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(iml_95),
      mean_y = mean(post_year),
      n = n(),
      .groups = 'drop'
    ) %>%
    filter(n >= 5)
  
  # Create RDD plot
  title_prefix <- ifelse(is.null(region_num), "National", paste("Region", region_num))
  
  rdd_plot <- ggplot(bin_means, aes(x = mean_x, y = mean_y)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_vline(xintercept = optimal_cutoff, 
               linetype = "dashed", 
               color = ifelse(optimal_result$pv[1] < 0.05, "red", "gray")) +
    geom_smooth(data = subset(bin_means, mean_x < optimal_cutoff),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(bin_means, mean_x >= optimal_cutoff),
                method = "lm", se = TRUE, color = "darkblue") +
    labs(
      title = paste(title_prefix, "- Entry Timing RD Plot"),
      subtitle = paste("Threshold Year:", threshold_year,
                       "| Cutoff:", round(optimal_cutoff, 2), 
                       "| Effect:", round(optimal_result$coef[1], 3)),
      x = "Marginality Index (1995)",
      y = paste("Post", threshold_year, "Entry (0/1)")
    ) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          plot.subtitle = element_text(size = 8))
  
  return(list(
    region = ifelse(is.null(region_num), "National", region_num),
    threshold_year = threshold_year,
    optimal_cutoff = optimal_cutoff,
    n_obs = nrow(data),
    rd_estimate = optimal_result$coef[1],
    std_error = optimal_result$se[1],
    p_value = optimal_result$pv[1],
    significant = optimal_result$pv[1] < 0.05,
    n_left = optimal_result$N_h[1],
    n_right = optimal_result$N_h[2],
    rd_result = optimal_result,
    rdd_plot = rdd_plot
  ))
}

# Get region-specific cutoffs
get_region_cutoffs <- function(data, threshold_year) {
  regions <- unique(data$region41)
  cutoffs <- data.frame()
  
  for (region in regions) {
    region_data <- data %>% filter(region41 == region)
    if (nrow(region_data) < 30) next
    
    region_results <- analyze_region_rdd(region_data, region_num = region, 
                                         threshold_year = threshold_year)
    
    if (!is.null(region_results)) {
      cutoffs <- bind_rows(cutoffs, 
                           data.frame(region = region, 
                                      cutoff = region_results$optimal_cutoff))
    }
  }
  
  return(cutoffs)
}

# Normalized national RDD
normalized_national_rdd <- function(data, threshold_year) {
  regional_cutoffs <- get_region_cutoffs(data, threshold_year)
  
  # Fix the join - regional_cutoffs has "region" column, data has "region41"
  data <- data %>%
    mutate(post_year = ifelse(entry_year_30 >= threshold_year, 1, 0)) %>%
    left_join(regional_cutoffs, by = "region41")
    filter(!is.na(cutoff)) %>%
    mutate(normalized_iml = iml_95 - cutoff)
  
  if (nrow(data) < 30) return(NULL)
  
  rd_result <- rdrobust(
    y = data$post_year,
    x = data$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  # Create plot
  bin_data <- data %>%
    mutate(bin = cut_width(normalized_iml, width = 0.1, center = 0))
  
  bin_means <- bin_data %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(normalized_iml),
      mean_y = mean(post_year),
      n = n(),
      .groups = 'drop'
    ) %>%
    filter(n >= 5)
  
  norm_plot <- ggplot(bin_means, aes(x = mean_x, y = mean_y)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(data = subset(bin_means, mean_x < 0),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(bin_means, mean_x >= 0),
                method = "lm", se = TRUE, color = "darkblue") +
    labs(
      title = "Normalized National RDD: Program Entry Timing",
      subtitle = paste0("Threshold Year: ", threshold_year),
      x = "Normalized Marginality Index (1995)",
      y = paste0("Post-", threshold_year, " Entry (0/1)")
    ) +
    theme_minimal()
  
  return(list(
    threshold_year = threshold_year,
    rd_result = rd_result,
    normalized_data = data,
    n_obs = nrow(data),
    rd_estimate = rd_result$coef[1],
    std_error = rd_result$se[1],
    p_value = rd_result$pv[1],
    significant = rd_result$pv[1] < 0.05,
    plot = norm_plot
  ))
}

# ============================================================================
# MAIN MULTI-YEAR ANALYSIS
# ============================================================================

run_multi_year_analysis <- function(data, years_range = 1998:2006) {
  
  # Initialize storage
  regional_results <- data.frame()
  national_results <- data.frame()
  regional_plots <- list()
  national_plots <- list()
  
  all_regions <- sort(unique(data$region41))
  
  # Analyze each year
  for (year in years_range) {
    cat("\nAnalyzing threshold year:", year, "\n")
    
    # Regional analysis
    for (region in all_regions) {
      cat("  Processing Region", region, "...\n")
      
      region_results <- analyze_region_rdd(data, region_num = region, 
                                           threshold_year = year)
      
      if (!is.null(region_results)) {
        regional_results <- rbind(regional_results,
                                  data.frame(
                                    threshold_year = year,
                                    region = region,
                                    optimal_cutoff = region_results$optimal_cutoff,
                                    n_obs = region_results$n_obs,
                                    rd_estimate = region_results$rd_estimate,
                                    std_error = region_results$std_error,
                                    p_value = region_results$p_value,
                                    significant = region_results$significant
                                  ))
        
        regional_plots[[paste(year, region, sep = "_")]] <- region_results$rdd_plot
      }
    }
    
    # National analysis
    national_result <- normalized_national_rdd(data, threshold_year = year)
    
    if (!is.null(national_result)) {
      national_results <- rbind(national_results,
                                data.frame(
                                  threshold_year = year,
                                  n_obs = national_result$n_obs,
                                  rd_estimate = national_result$rd_estimate,
                                  std_error = national_result$std_error,
                                  p_value = national_result$p_value,
                                  significant = national_result$significant
                                ))
      
      national_plots[[as.character(year)]] <- national_result$plot
    }
  }
  
  # Save regional plots by year
  for (year in years_range) {
    year_plots <- regional_plots[grep(paste0("^", year, "_"), names(regional_plots))]
    
    if (length(year_plots) > 0) {
      n_plots <- length(year_plots)
      n_cols <- min(3, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      grid_plot <- do.call(grid.arrange, 
                           c(year_plots, 
                             list(ncol = n_cols, nrow = n_rows,
                                  top = paste("Regional RDD Analysis - Threshold Year:", year))))
      
      ggsave(paste0("multi_year_analysis/regional_plots/year_", year, "_regions.pdf"), 
             grid_plot, width = 15, height = 12)
    }
  }
  
  # Create national plots grid
  if (length(national_plots) > 0) {
    n_plots <- length(national_plots)
    n_cols <- min(3, n_plots)
    n_rows <- ceiling(n_plots / n_cols)
    
    national_grid <- do.call(grid.arrange, 
                             c(national_plots, 
                               list(ncol = n_cols, nrow = n_rows,
                                    top = "Normalized National RDD Analysis by Threshold Year")))
    
    ggsave("multi_year_analysis/national_rdd_all_years.pdf", 
           national_grid, width = 15, height = 12)
  }
  
  # Create results tables with significance stars
  regional_table_wide <- regional_results %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1 ~ "*",
        TRUE ~ ""
      ),
      result = paste0(round(rd_estimate, 3), significance)
    ) %>%
    select(threshold_year, region, result) %>%
    pivot_wider(names_from = threshold_year, values_from = result)
  
  national_table <- national_results %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.1 ~ "*",
        TRUE ~ ""
      ),
      result = paste0(round(rd_estimate, 3), significance),
      threshold_year = as.character(threshold_year)
    ) %>%
    select(threshold_year, result, n_obs) %>%
    arrange(threshold_year)
  
  # National trend visualization
  national_trend <- ggplot(national_results, aes(x = factor(threshold_year), y = rd_estimate)) +
    geom_point(aes(size = 1/std_error, color = p_value < 0.05)) +
    geom_errorbar(aes(ymin = rd_estimate - 1.96 * std_error, 
                      ymax = rd_estimate + 1.96 * std_error, 
                      color = p_value < 0.05), 
                  width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("gray", "blue"), name = "Significant (p<0.05)") +
    labs(
      title = "National RDD Estimates by Threshold Year",
      subtitle = "Using Normalized Regional Cutoffs",
      x = "Threshold Year",
      y = "RD Estimate with 95% CI",
      size = "Precision (1/SE)"
    ) +
    theme_minimal()
  
  ggsave("multi_year_analysis/national_trend.pdf", national_trend, width = 10, height = 6)
  
  # LaTeX tables
  stargazer(
    regional_table_wide, 
    type = "latex", 
    summary = FALSE, 
    rownames = FALSE,
    title = "Regional RDD Results by Entry Year Threshold", 
    digits = 3, 
    align = TRUE, 
    out = "multi_year_analysis/regional_results_table.tex"
  )
  
  stargazer(
    national_table, 
    type = "latex", 
    summary = FALSE, 
    rownames = FALSE,
    title = "Normalized National RDD Results by Entry Year Threshold", 
    digits = 3, 
    align = TRUE, 
    out = "multi_year_analysis/national_results_table.tex"
  )
  
  return(list(
    regional_results = regional_results,
    national_results = national_results,
    regional_table = regional_table_wide,
    national_table = national_table,
    regional_plots = regional_plots,
    national_plots = national_plots,
    national_trend = national_trend
  ))
}

# ============================================================================
# RUN MAIN ANALYSIS
# ============================================================================

cat("Running multi-year analysis for years 1998 to 2006...\n")
results <- run_multi_year_analysis(merged_entry_iml_rg, years_range = 1998:2006)

# Print summary
cat("\n========================================\n")
cat("Summary of National Results:\n")
cat("========================================\n")
national_summary <- results$national_results %>%
  select(threshold_year, rd_estimate, std_error, p_value, significant) %>%
  mutate_if(is.numeric, round, 3) %>%
  arrange(threshold_year)

print(national_summary)

# Save results
save(results, file = "multi_year_analysis/rdd_results_all_years.RData")

# ============================================================================
# CREATE REGIONAL TIME SERIES PLOTS (matching original style)
# ============================================================================

cat("\nCreating regional time series plots...\n")

regional_time_series <- results$regional_results %>%
  select(threshold_year, region, rd_estimate, significant) %>%
  mutate(
    year = as.numeric(threshold_year),
    region = as.factor(region)
  )

# Create small multiple plots for each region
create_region_time_plots <- function() {
  all_regions <- sort(unique(regional_time_series$region))
  plot_list <- list()
  
  # Vertical line for 2001-2002 administrative change
  vertical_line <- geom_vline(xintercept = 2001.5, color = "red", 
                              linetype = "dashed", alpha = 0.7)
  
  for (r in all_regions) {
    region_data <- regional_time_series %>% 
      filter(region == r)
    
    if (nrow(region_data) < 2) next
    
    p <- ggplot(region_data, aes(x = year, y = rd_estimate)) +
      vertical_line +
      geom_line(color = "blue", alpha = 0.7) +
      geom_point(aes(shape = significant), size = 2.5) +
      scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8)) +
      scale_x_continuous(breaks = seq(2000, 2006, by = 2)) +
      labs(title = r) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 9, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(size = 6),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
      )
    
    plot_list[[as.character(r)]] <- p
  }
  
  return(plot_list)
}

region_plots <- create_region_time_plots()

# Arrange in grid
if (length(region_plots) > 0) {
  n_plots <- length(region_plots)
  n_cols <- 6
  n_rows <- ceiling(n_plots / n_cols)
  
  # Create legend
  legend_data <- data.frame(
    year = 2002,
    rd_estimate = 0,
    significant = c(TRUE, FALSE)
  )
  
  legend_plot <- ggplot(legend_data, aes(x = year, y = rd_estimate, shape = significant)) +
    geom_point() +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 8),
                       labels = c("FALSE" = "Not significant", "TRUE" = "Significant at 5%")) +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  legend <- cowplot::get_legend(legend_plot)
  
  # Title elements
  main_title <- textGrob("RD Estimates Over Time by Region", 
                         gp = gpar(fontsize = 14, fontface = "bold"))
  subtitle <- textGrob("Red dashed line indicates 2001-2002 administrative change", 
                       gp = gpar(fontsize = 10))
  
  # Create final grid
  region_grid <- arrangeGrob(
    grobs = region_plots,
    ncol = n_cols,
    nrow = n_rows,
    top = arrangeGrob(main_title, subtitle, ncol = 1, heights = c(1, 0.6)),
    bottom = legend
  )
  
  x_axis_label <- textGrob("Year", gp = gpar(fontsize = 10))
  y_axis_label <- textGrob("RD Estimate", gp = gpar(fontsize = 10), rot = 90)
  
  final_grid <- grid.arrange(
    arrangeGrob(region_grid, bottom = x_axis_label, left = y_axis_label)
  )
  
  ggsave("multi_year_analysis/rd_estimates_over_time_by_region.pdf", 
         final_grid, width = 14, height = 10)
  ggsave("multi_year_analysis/rd_estimates_over_time_by_region.png", 
         final_grid, width = 14, height = 10, dpi = 300)
}

# ============================================================================
# CREATE HEATMAP VISUALIZATION
# ============================================================================

regional_heatmap_data <- results$regional_results %>%
  select(threshold_year, region, rd_estimate, significant) %>%
  mutate(
    threshold_year = as.factor(threshold_year),
    region = as.factor(region),
    rd_estimate = round(rd_estimate, 3),
    display_value = ifelse(significant, paste0(rd_estimate, "*"), as.character(rd_estimate))
  )

heatmap_plot <- ggplot(regional_heatmap_data, 
                       aes(x = threshold_year, y = region, fill = rd_estimate)) +
  geom_tile() +
  geom_text(aes(label = display_value), size = 2.5) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, name = "RD Estimate") +
  labs(
    title = "RD Estimates by Region and Threshold Year",
    subtitle = "* indicates statistical significance at 5% level",
    x = "Entry Year Threshold",
    y = "Region"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("multi_year_analysis/regional_heatmap.pdf", heatmap_plot, width = 12, height = 10)

# ============================================================================
# PEAK EFFECT YEAR ANALYSIS
# ============================================================================

max_effect_year <- results$national_results %>%
  filter(abs(rd_estimate) == max(abs(rd_estimate))) %>%
  pull(threshold_year)

cat("\nYear with largest effect:", max_effect_year, "\n")

if (length(max_effect_year) > 0) {
  max_year_result <- normalized_national_rdd(merged_entry_iml_rg, 
                                             threshold_year = max_effect_year)
  
  if (!is.null(max_year_result)) {
    before_after_plot <- ggplot(max_year_result$normalized_data) +
      geom_density(aes(x = normalized_iml, fill = factor(post_year)), 
                   alpha = 0.5, position = "identity") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
      scale_fill_manual(values = c("turquoise", "coral"),
                        name = "Entry Period",
                        labels = c(paste0("Pre-", max_effect_year), 
                                   paste0("Post-", max_effect_year))) +
      labs(
        title = paste0("Distribution of Normalized IML by Entry Period",
                       " (Threshold Year: ", max_effect_year, ")"),
        subtitle = paste0("RD Estimate: ", round(max_year_result$rd_estimate, 3), 
                          " | p-value: ", round(max_year_result$p_value, 3)),
        x = "Normalized Marginality Index (1995)",
        y = "Density"
      ) +
      theme_minimal()
    
    ggsave("multi_year_analysis/peak_effect_year_distribution.pdf", 
           before_after_plot, width = 10, height = 6)
  }
}

cat("\n========================================\n")
cat("TEMPORAL ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("\nAll outputs saved in multi_year_analysis/\n")

# Override the functions with debug versions
get_region_cutoffs <- function(data, threshold_year) {
  regions <- unique(data$region41)
  cutoffs <- data.frame()
  
  for (region in regions) {
    region_data <- data %>% filter(region41 == region)
    if (nrow(region_data) < 30) next
    
    region_results <- analyze_region_rdd(region_data, region_num = region, 
                                        threshold_year = threshold_year)
    
    if (!is.null(region_results)) {
      cutoffs <- bind_rows(cutoffs, 
                          data.frame(region = region,  # Keep as "region" for now
                                   cutoff = region_results$optimal_cutoff))
    }
  }
  
  print(paste("Cutoffs columns:", paste(names(cutoffs), collapse=", ")))  # DEBUG
  return(cutoffs)
}

normalized_national_rdd <- function(data, threshold_year) {
  regional_cutoffs <- get_region_cutoffs(data, threshold_year)
  
  print(paste("Data has columns:", paste(names(data)[1:10], collapse=", "), "..."))  # DEBUG
  print(paste("Cutoffs has columns:", paste(names(regional_cutoffs), collapse=", ")))  # DEBUG
  
  # Rename column to match
  if("region" %in% names(regional_cutoffs)) {
    regional_cutoffs <- regional_cutoffs %>% rename(region41 = region)
  }
  
  data <- data %>%
    mutate(post_year = ifelse(entry_year_30 >= threshold_year, 1, 0)) %>%
    left_join(regional_cutoffs, by = "region41") %>%
    filter(!is.na(cutoff)) %>%
    mutate(normalized_iml = iml_95 - cutoff)
  
  if (nrow(data) < 30) return(NULL)
  
  rd_result <- rdrobust(
    y = data$post_year,
    x = data$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  # ... rest of function stays the same
  
  bin_data <- data %>%
    mutate(bin = cut_width(normalized_iml, width = 0.1, center = 0))
  
  bin_means <- bin_data %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(normalized_iml),
      mean_y = mean(post_year),
      n = n(),
      .groups = 'drop'
    ) %>%
    filter(n >= 5)
  
  norm_plot <- ggplot(bin_means, aes(x = mean_x, y = mean_y)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(data = subset(bin_means, mean_x < 0),
                method = "lm", se = TRUE, color = "darkblue") +
    geom_smooth(data = subset(bin_means, mean_x >= 0),
                method = "lm", se = TRUE, color = "darkblue") +
    labs(
      title = "Normalized National RDD: Program Entry Timing",
      subtitle = paste0("Threshold Year: ", threshold_year),
      x = "Normalized Marginality Index (1995)",
      y = paste0("Post-", threshold_year, " Entry (0/1)")
    ) +
    theme_minimal()
  
  return(list(
    threshold_year = threshold_year,
    rd_result = rd_result,
    normalized_data = data,
    n_obs = nrow(data),
    rd_estimate = rd_result$coef[1],
    std_error = rd_result$se[1],
    p_value = rd_result$pv[1],
    significant = rd_result$pv[1] < 0.05,
    plot = norm_plot
  ))
}