# 02_targeting_analysis.R
# Analysis of PUN_INC distributions by TIPO_CAP and state-level targeting patterns

# Load required libraries
library(tidyverse)
library(gridExtra)
library(scales)
library(viridis)
library(patchwork)

# Ensure scales functions are available
comma <- scales::comma

# Define color scheme
color_scheme <- c(
  "Rural Census" = "cornflowerblue",
  "Urban Census" = "steelblue",
  "Ventanilla" = "coral"
)

# Load prepared data
prepared_data <- readRDS("prepared_encaseh_data.rds")

# ========================================
# PART 1: PUN_INC DISTRIBUTIONS BY TIPO_CAP
# ========================================

analyze_pun_inc_by_tipo_cap <- function(data, year) {
  # Filter for valid TIPO_CAP and PUN_INC values
  df <- data %>%
    filter(!is.na(targeting_method) & !is.na(PUN_INC)) %>%
    filter(PUN_INC < 100)  # Focus on normal scoring system for cleaner visualization
  
  # Calculate summary statistics
  summary_stats <- df %>%
    group_by(targeting_method) %>%
    summarise(
      n = n(),
      mean_score = mean(PUN_INC, na.rm = TRUE),
      median_score = median(PUN_INC, na.rm = TRUE),
      sd_score = sd(PUN_INC, na.rm = TRUE),
      q25 = quantile(PUN_INC, 0.25, na.rm = TRUE),
      q75 = quantile(PUN_INC, 0.75, na.rm = TRUE),
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100
    ) %>%
    mutate(across(where(is.numeric), ~round(., 2)))
  
  # Create density plot
  p1 <- ggplot(df, aes(x = PUN_INC, fill = targeting_method, color = targeting_method)) +
    geom_density(alpha = 0.4, size = 1) +
    facet_wrap(~targeting_method, nrow = 3, scales = "free_y") +
    scale_fill_manual(values = color_scheme) +
    scale_color_manual(values = color_scheme) +
    labs(
      title = sprintf("PUN_INC Distribution by Targeting Method - %d", year),
      x = "Poverty Score (PUN_INC)",
      y = "Density",
      subtitle = "Excluding high-value outliers (PUN_INC > 100)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA)
    )
  
  # Create boxplot comparison
  p2 <- ggplot(df, aes(x = targeting_method, y = PUN_INC, fill = targeting_method)) +
    geom_boxplot(alpha = 0.6, outlier.alpha = 0.3) +
    geom_violin(alpha = 0.3) +
    scale_fill_manual(values = color_scheme) +
    coord_flip() +
    labs(
      title = "PUN_INC Distribution Comparison",
      x = "",
      y = "Poverty Score (PUN_INC)",
      subtitle = sprintf("n = %s observations", format(nrow(df), big.mark = ","))
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Create acceptance rate by score bins
  df_bins <- df %>%
    mutate(score_bin = cut(PUN_INC, 
                           breaks = quantile(PUN_INC, probs = seq(0, 1, 0.1), na.rm = TRUE),
                           include.lowest = TRUE,
                           labels = paste0("D", 1:10))) %>%
    group_by(targeting_method, score_bin) %>%
    summarise(
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
      n = n(),
      .groups = "drop"
    )
  
  p3 <- ggplot(df_bins, aes(x = score_bin, y = acceptance_rate, 
                            color = targeting_method, group = targeting_method)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    scale_color_manual(values = color_scheme) +
    labs(
      title = "Acceptance Rate by Poverty Score Decile",
      x = "Poverty Score Decile (D1 = poorest)",
      y = "Acceptance Rate (%)",
      color = "Targeting Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Combine plots using patchwork
  combined_plot <- p1 / p2 / p3
  
  combined_plot <- combined_plot + 
    plot_annotation(
      title = sprintf("Poverty Score Analysis by Targeting Method - Year %d", year),
      theme = theme(plot.title = element_text(size = 16, face = "bold"))
    )
  
  return(list(
    plot = combined_plot,
    summary = summary_stats
  ))
}

# ========================================
# PART 2: STATE-LEVEL TARGETING PATTERNS
# ========================================

analyze_state_targeting <- function(data, year) {
  # Get state names mapping
  state_names <- c(
    "1" = "Aguascalientes", "2" = "Baja California", "3" = "Baja California Sur",
    "4" = "Campeche", "5" = "Coahuila", "6" = "Colima", "7" = "Chiapas",
    "8" = "Chihuahua", "9" = "Ciudad de México", "10" = "Durango",
    "11" = "Guanajuato", "12" = "Guerrero", "13" = "Hidalgo", "14" = "Jalisco",
    "15" = "México", "16" = "Michoacán", "17" = "Morelos", "18" = "Nayarit",
    "19" = "Nuevo León", "20" = "Oaxaca", "21" = "Puebla", "22" = "Querétaro",
    "23" = "Quintana Roo", "24" = "San Luis Potosí", "25" = "Sinaloa",
    "26" = "Sonora", "27" = "Tabasco", "28" = "Tamaulipas", "29" = "Tlaxcala",
    "30" = "Veracruz", "31" = "Yucatán", "32" = "Zacatecas"
  )
  
  # Filter for valid TIPO_CAP
  df <- data %>%
    filter(!is.na(targeting_method) & !is.na(ENT_FED)) %>%
    mutate(state_name = state_names[as.character(ENT_FED)])
  
  # 1. Number of localities using each method by state
  localities_by_method <- df %>%
    group_by(ENT_FED, state_name, LOCALI, targeting_method) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ENT_FED, state_name, targeting_method) %>%
    summarise(n_localities = n_distinct(LOCALI), .groups = "drop") %>%
    pivot_wider(names_from = targeting_method, values_from = n_localities, values_fill = 0) %>%
    pivot_longer(cols = -c(ENT_FED, state_name), names_to = "targeting_method", values_to = "n_localities")
  
  # 2. Total households by method and state
  households_by_method <- df %>%
    group_by(ENT_FED, state_name, targeting_method) %>%
    summarise(
      n_households = n(),
      n_accepted = sum(accepted, na.rm = TRUE),
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Create visualization for localities
  p1 <- ggplot(localities_by_method %>% filter(n_localities > 0), 
               aes(x = reorder(state_name, ENT_FED), y = n_localities, fill = targeting_method)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = color_scheme) +
    coord_flip() +
    labs(
      title = sprintf("Number of Localities by Targeting Method - %d", year),
      x = "",
      y = "Number of Localities",
      fill = "Targeting Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  # Create visualization for households
  p2 <- ggplot(households_by_method, 
               aes(x = reorder(state_name, ENT_FED), y = n_households, fill = targeting_method)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = color_scheme) +
    scale_y_continuous(labels = scales::comma) +
    coord_flip() +
    labs(
      title = sprintf("Number of Households by Targeting Method - %d", year),
      x = "",
      y = "Number of Households",
      fill = "Targeting Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  # Create percentage composition plot
  household_pct <- households_by_method %>%
    group_by(ENT_FED, state_name) %>%
    mutate(
      total_households = sum(n_households),
      pct_households = (n_households / total_households) * 100
    ) %>%
    ungroup()
  
  p3 <- ggplot(household_pct, 
               aes(x = reorder(state_name, ENT_FED), y = pct_households, fill = targeting_method)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = color_scheme) +
    coord_flip() +
    labs(
      title = "Percentage Distribution of Targeting Methods by State",
      x = "",
      y = "Percentage of Households (%)",
      fill = "Targeting Method"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 8)
    )
  
  # Create acceptance rate heatmap
  p4 <- ggplot(households_by_method, 
               aes(x = targeting_method, y = reorder(state_name, -ENT_FED), fill = acceptance_rate)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = sprintf("%.1f%%", acceptance_rate)), size = 3) +
    scale_fill_gradient2(low = "lightblue", mid = "white", high = "darkred", 
                         midpoint = 50, limits = c(0, 100)) +
    labs(
      title = "Acceptance Rates by State and Targeting Method",
      x = "Targeting Method",
      y = "",
      fill = "Acceptance\nRate (%)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
  
  # Summary statistics
  summary_stats <- households_by_method %>%
    group_by(targeting_method) %>%
    summarise(
      total_states = n_distinct(ENT_FED),
      total_households = sum(n_households),
      total_accepted = sum(n_accepted),
      avg_acceptance_rate = weighted.mean(acceptance_rate, n_households),
      .groups = "drop"
    )
  
  return(list(
    localities_plot = p1,
    households_plot = p2,
    percentage_plot = p3,
    heatmap_plot = p4,
    summary = summary_stats
  ))
}

# ========================================
# PART 3: TEMPORAL ANALYSIS
# ========================================

create_temporal_comparison <- function(prepared_data) {
  # Combine data from years with TIPO_CAP
  years_with_tipo <- list(
    "2001" = prepared_data$y2001,
    "2002" = prepared_data$y2002,
    "2004" = prepared_data$y2004
  )
  
  temporal_summary <- map_df(names(years_with_tipo), function(yr) {
    data <- years_with_tipo[[yr]]
    
    data %>%
      filter(!is.na(targeting_method)) %>%
      group_by(targeting_method) %>%
      summarise(
        n_households = n(),
        n_accepted = sum(accepted, na.rm = TRUE),
        acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
        mean_score = mean(PUN_INC[PUN_INC < 100], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(year = as.numeric(yr))
  })
  
  # Plot temporal evolution
  p1 <- ggplot(temporal_summary, aes(x = year, y = n_households, 
                                     color = targeting_method, group = targeting_method)) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    scale_color_manual(values = color_scheme) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title = "Evolution of Targeting Methods Over Time",
      x = "Year",
      y = "Number of Households",
      color = "Targeting Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  p2 <- ggplot(temporal_summary, aes(x = year, y = acceptance_rate, 
                                     color = targeting_method, group = targeting_method)) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    scale_color_manual(values = color_scheme) +
    labs(
      title = "Acceptance Rates by Targeting Method Over Time",
      x = "Year",
      y = "Acceptance Rate (%)",
      color = "Targeting Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  combined_temporal <- p1 / p2
  
  return(list(
    plot = combined_temporal,
    summary = temporal_summary
  ))
}

# ========================================
# RUN ANALYSIS
# ========================================

# Create output directory if it doesn't exist
dir.create("output", showWarnings = FALSE)

# Analyze each year with TIPO_CAP
cat("Analyzing PUN_INC distributions by TIPO_CAP...\n")

# 2001 Analysis
pun_inc_2001 <- analyze_pun_inc_by_tipo_cap(prepared_data$y2001, 2001)
state_2001 <- analyze_state_targeting(prepared_data$y2001, 2001)

# 2002 Analysis
pun_inc_2002 <- analyze_pun_inc_by_tipo_cap(prepared_data$y2002, 2002)
state_2002 <- analyze_state_targeting(prepared_data$y2002, 2002)

# 2004 Analysis
pun_inc_2004 <- analyze_pun_inc_by_tipo_cap(prepared_data$y2004, 2004)
state_2004 <- analyze_state_targeting(prepared_data$y2004, 2004)

# Temporal comparison
temporal_analysis <- create_temporal_comparison(prepared_data)

# Save PUN_INC distribution plots as PDF
ggsave("output/pun_inc_distribution_2001.pdf", pun_inc_2001$plot, 
       width = 12, height = 16)
ggsave("output/pun_inc_distribution_2002.pdf", pun_inc_2002$plot, 
       width = 12, height = 16)
ggsave("output/pun_inc_distribution_2004.pdf", pun_inc_2004$plot, 
       width = 12, height = 16)

# Save state-level plots for ALL years as PDF
# 2001
ggsave("output/state_localities_2001.pdf", state_2001$localities_plot, 
       width = 10, height = 12)
ggsave("output/state_households_2001.pdf", state_2001$households_plot, 
       width = 10, height = 12)
ggsave("output/state_percentage_2001.pdf", state_2001$percentage_plot, 
       width = 10, height = 12)
ggsave("output/state_heatmap_2001.pdf", state_2001$heatmap_plot, 
       width = 10, height = 12)

# 2002
ggsave("output/state_localities_2002.pdf", state_2002$localities_plot, 
       width = 10, height = 12)
ggsave("output/state_households_2002.pdf", state_2002$households_plot, 
       width = 10, height = 12)
ggsave("output/state_percentage_2002.pdf", state_2002$percentage_plot, 
       width = 10, height = 12)
ggsave("output/state_heatmap_2002.pdf", state_2002$heatmap_plot, 
       width = 10, height = 12)

# 2004
ggsave("output/state_localities_2004.pdf", state_2004$localities_plot, 
       width = 10, height = 12)
ggsave("output/state_households_2004.pdf", state_2004$households_plot, 
       width = 10, height = 12)
ggsave("output/state_percentage_2004.pdf", state_2004$percentage_plot, 
       width = 10, height = 12)
ggsave("output/state_heatmap_2004.pdf", state_2004$heatmap_plot, 
       width = 10, height = 12)

# Save temporal evolution as PDF
ggsave("output/temporal_evolution.pdf", temporal_analysis$plot, 
       width = 10, height = 10)

# Print summaries
cat("\n=== PUN_INC Distribution Summary 2002 ===\n")
print(pun_inc_2002$summary)

cat("\n=== State-Level Targeting Summary 2002 ===\n")
print(state_2002$summary)

cat("\n=== Temporal Evolution Summary ===\n")
print(temporal_analysis$summary)

cat("\n=== PUN_INC Distribution Summary 2004 ===\n")
print(pun_inc_2004$summary)

cat("\n=== State-Level Targeting Summary 2004 ===\n")
print(state_2004$summary)

# ========================================
# ADDITIONAL ANALYSIS: Focus on Ventanilla
# ========================================

analyze_ventanilla_expansion <- function(prepared_data) {
  # Identify states with significant Ventanilla presence
  ventanilla_states <- map_df(c("2001", "2002", "2004"), function(yr) {
    data <- switch(yr,
                   "2001" = prepared_data$y2001,
                   "2002" = prepared_data$y2002,
                   "2004" = prepared_data$y2004
    )
    
    data %>%
      filter(targeting_method == "Ventanilla") %>%
      group_by(ENT_FED) %>%
      summarise(
        n_ventanilla = n(),
        n_localities = n_distinct(LOCALI),
        mean_score = mean(PUN_INC[PUN_INC < 100], na.rm = TRUE),
        acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
        .groups = "drop"
      ) %>%
      mutate(year = as.numeric(yr))
  })
  
  # Top states using Ventanilla
  top_ventanilla_states <- ventanilla_states %>%
    group_by(ENT_FED) %>%
    summarise(total_ventanilla = sum(n_ventanilla), .groups = "drop") %>%
    arrange(desc(total_ventanilla)) %>%
    slice_head(n = 10)
  
  # Create focused visualization
  p <- ventanilla_states %>%
    filter(ENT_FED %in% top_ventanilla_states$ENT_FED) %>%
    ggplot(aes(x = factor(year), y = n_ventanilla, group = factor(ENT_FED))) +
    geom_line(aes(color = factor(ENT_FED)), size = 1.2) +
    geom_point(aes(color = factor(ENT_FED)), size = 3) +
    scale_y_continuous(labels = comma) +
    labs(
      title = "Ventanilla (Self-Selection) Expansion in Top 10 States",
      subtitle = "States with highest adoption of self-selection mechanism",
      x = "Year",
      y = "Number of Households using Ventanilla",
      color = "State Code"
    ) +
    theme_minimal() +
    facet_wrap(~ENT_FED, scales = "free_y", ncol = 2)
  
  return(list(plot = p, summary = top_ventanilla_states))
}

ventanilla_analysis <- analyze_ventanilla_expansion(prepared_data)
ggsave("output/ventanilla_expansion.pdf", ventanilla_analysis$plot, 
       width = 12, height = 10)

cat("\n=== Top States Using Ventanilla ===\n")
print(ventanilla_analysis$summary)

cat("\nAnalysis complete! Check the 'output' folder for PDF visualizations.\n")
