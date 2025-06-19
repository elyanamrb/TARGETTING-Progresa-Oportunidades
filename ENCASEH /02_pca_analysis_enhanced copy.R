# 02_pca_analysis_enhanced.R
# Enhanced PCA analysis that separates the two scoring systems

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(kableExtra)

# Load prepared data
prepared_data <- readRDS("prepared_encaseh_data.rds")

# Define asset variables
asset_vars <- c("LICUADORA", "REFRI", "EST_GAS", "CAL_AGUA", "RADIO", 
                "TOCADISC", "TELEVIS", "VIDEOCAS", "LAVADORA", "VENTILA", 
                "AUTO_PROP", "CAMION_PRO", "TRACTOR")

# Create directories for outputs
dir.create("tables", showWarnings = FALSE)

# Enhanced function to identify and analyze scoring systems
analyze_scoring_systems <- function(data, year_label) {
  cat(sprintf("\n=== Analyzing Scoring Systems for Year %s ===\n", year_label))
  
  # Identify scoring systems based on PUN_INC distribution
  data <- data %>%
    mutate(
      scoring_system_detailed = case_when(
        PUN_INC > 1000 ~ "System 2 (Very High >1000)",
        PUN_INC > 100 ~ "System 2 (High 100-1000)",
        PUN_INC > 10 ~ "System 1 (Normal 10-100)",
        PUN_INC >= 0 ~ "System 1 (Low 0-10)",
        TRUE ~ "Unknown"
      ),
      scoring_system_binary = ifelse(PUN_INC > 100, "System 2", "System 1")
    )
  
  # Summary statistics by scoring system
  system_summary <- data %>%
    group_by(scoring_system_binary) %>%
    summarise(
      n = n(),
      pct = n() / nrow(data) * 100,
      mean_pun_inc = mean(PUN_INC, na.rm = TRUE),
      median_pun_inc = median(PUN_INC, na.rm = TRUE),
      sd_pun_inc = sd(PUN_INC, na.rm = TRUE),
      min_pun_inc = min(PUN_INC, na.rm = TRUE),
      max_pun_inc = max(PUN_INC, na.rm = TRUE),
      acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
      .groups = "drop"
    )
  
  # Regional distribution of scoring systems
  regional_systems <- data %>%
    filter(!is.na(REGION)) %>%
    group_by(REGION, scoring_system_binary) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = scoring_system_binary, values_from = n, values_fill = 0) %>%
    mutate(
      total = rowSums(select(., -REGION), na.rm = TRUE),
      pct_system2 = ifelse("System 2" %in% names(.), `System 2` / total * 100, 0)
    ) %>%
    arrange(desc(pct_system2))
  
  # Identify regions using predominantly System 2
  system2_regions <- regional_systems %>%
    filter(pct_system2 > 50) %>%
    pull(REGION)
  
  # Create visualizations
  # 1. Distribution plot by scoring system (handle negative values)
  data_for_plot <- data %>%
    filter(PUN_INC > 0)  # Remove non-positive values for log scale
  
  dist_plot <- ggplot(data_for_plot, aes(x = PUN_INC, fill = scoring_system_binary)) +
    geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
    scale_x_log10(labels = scales::comma) +
    facet_wrap(~scoring_system_binary, scales = "free") +
    labs(title = sprintf("PUN_INC Distribution by Scoring System - %s", year_label),
         subtitle = sprintf("System 1: %.1f%% of households, System 2: %.1f%% of households (excluding %d non-positive values)",
                            sum(data$scoring_system_binary == "System 1", na.rm = TRUE) / nrow(data) * 100,
                            sum(data$scoring_system_binary == "System 2", na.rm = TRUE) / nrow(data) * 100,
                            sum(data$PUN_INC <= 0, na.rm = TRUE)),
         x = "PUN_INC (log scale)",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 2. Regional heatmap
  if (length(unique(data$REGION[!is.na(data$REGION)])) > 1) {
    regional_plot <- regional_systems %>%
      filter(total > 10) %>%  # Only show regions with sufficient data
      ggplot(aes(x = reorder(factor(REGION), pct_system2), y = pct_system2)) +
      geom_bar(stat = "identity", fill = "coral") +
      geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
      coord_flip() +
      labs(title = sprintf("Percentage Using System 2 by Region - %s", year_label),
           x = "Region",
           y = "% Using System 2 (High scores)") +
      theme_minimal()
  } else {
    regional_plot <- NULL
  }
  
  # Save plots
  ggsave(sprintf("descriptive_stats/scoring_systems_dist_%s.pdf", year_label), 
         dist_plot, width = 12, height = 6)
  
  if (!is.null(regional_plot)) {
    ggsave(sprintf("descriptive_stats/scoring_systems_regional_%s.pdf", year_label), 
           regional_plot, width = 10, height = 12)
  }
  
  return(list(
    data = data,
    system_summary = system_summary,
    regional_systems = regional_systems,
    system2_regions = system2_regions,
    plots = list(dist_plot = dist_plot, regional_plot = regional_plot)
  ))
}

# Enhanced PCA function that handles scoring systems separately
perform_pca_by_system <- function(data, year_label, scoring_analysis) {
  cat(sprintf("\n=== PCA Analysis by Scoring System for Year %s ===\n", year_label))
  
  # Use the data with scoring system assignments
  data <- scoring_analysis$data
  
  # Check which assets exist in this year's data
  existing_assets <- asset_vars[asset_vars %in% names(data)]
  
  # Initialize results storage
  pca_results_by_system <- list()
  
  # Perform PCA for each scoring system
  for (system in unique(data$scoring_system_binary)) {
    cat(sprintf("\n--- Analyzing %s ---\n", system))
    
    # Filter data for this system
    system_data <- data %>% filter(scoring_system_binary == system)
    
    # Prepare asset data
    assets_df <- system_data[, existing_assets]
    
    # Recode assets: 1 = has asset, 0 = doesn't have asset
    for(col in names(assets_df)) {
      assets_df[[col]] <- case_when(
        assets_df[[col]] == 1 ~ 1,
        assets_df[[col]] == 2 ~ 0,
        TRUE ~ NA_real_
      )
    }
    
    # Get complete cases
    assets_complete <- na.omit(assets_df)
    n_complete <- nrow(assets_complete)
    cat(sprintf("Complete cases: %d out of %d (%.1f%%)\n", 
                n_complete, nrow(system_data), n_complete/nrow(system_data)*100))
    
    # Skip if too few complete cases
    if (n_complete < 50) {
      cat("Too few complete cases for PCA analysis\n")
      pca_results_by_system[[system]] <- NULL
      next
    }
    
    # Run PCA
    pca_result <- prcomp(assets_complete, center = TRUE, scale. = TRUE)
    
    # Get PC1 scores
    pc1_scores <- pca_result$x[,1]
    
    # Match back to original data
    complete_indices <- as.numeric(rownames(assets_complete))
    
    # Create correlation dataset
    correlation_data <- data.frame(
      pun_inc = system_data$PUN_INC[complete_indices],
      pca_score = pc1_scores,
      pca_score_reversed = -pc1_scores  # Reverse so higher = more assets
    ) %>% na.omit()
    
    # Calculate correlations
    cor_original <- cor(correlation_data$pun_inc, correlation_data$pca_score)
    cor_reversed <- cor(correlation_data$pun_inc, correlation_data$pca_score_reversed)
    
    # Store results
    pca_results_by_system[[system]] <- list(
      system = system,
      n_complete = n_complete,
      pca_result = pca_result,
      correlation = cor_reversed,
      correlation_data = correlation_data,
      variance_explained = pca_result$sdev^2 / sum(pca_result$sdev^2)
    )
    
    cat(sprintf("Correlation (reversed PCA vs PUN_INC): %.3f\n", cor_reversed))
    cat(sprintf("Variance explained by PC1: %.1f%%\n", 
                pca_results_by_system[[system]]$variance_explained[1] * 100))
  }
  
  # Create comparison visualizations
  if (length(pca_results_by_system) > 1) {
    # Combined scatter plot
    scatter_data <- data.frame()
    for (system in names(pca_results_by_system)) {
      if (!is.null(pca_results_by_system[[system]])) {
        temp_data <- pca_results_by_system[[system]]$correlation_data %>%
          mutate(System = system)
        scatter_data <- rbind(scatter_data, temp_data)
      }
    }
    
    combined_scatter <- ggplot(scatter_data, aes(x = pca_score_reversed, y = pun_inc)) +
      geom_point(alpha = 0.3, aes(color = System), size = 0.8) +
      geom_smooth(method = "lm", aes(color = System), se = TRUE) +
      facet_wrap(~System, scales = "free") +
      labs(title = sprintf("PCA Score vs PUN_INC by Scoring System - %s", year_label),
           x = "PCA Score (Reversed - Higher = More Assets)",
           y = "PUN_INC") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Loadings comparison
    loadings_comparison <- data.frame()
    for (system in names(pca_results_by_system)) {
      if (!is.null(pca_results_by_system[[system]])) {
        pc1_loadings <- pca_results_by_system[[system]]$pca_result$rotation[,1]
        temp_loadings <- data.frame(
          Asset = names(pc1_loadings),
          Loading = pc1_loadings,
          System = system
        )
        loadings_comparison <- rbind(loadings_comparison, temp_loadings)
      }
    }
    
    loadings_plot <- ggplot(loadings_comparison, 
                            aes(x = reorder(Asset, Loading), y = Loading, fill = System)) +
      geom_bar(stat = "identity", position = "dodge") +
      coord_flip() +
      labs(title = sprintf("PCA Loadings Comparison by System - %s", year_label),
           x = "Asset",
           y = "Loading on PC1") +
      theme_minimal()
    
    # Save comparison plots
    ggsave(sprintf("pca_plots/pca_scatter_by_system_%s.pdf", year_label), 
           combined_scatter, width = 12, height = 6)
    ggsave(sprintf("pca_plots/pca_loadings_by_system_%s.pdf", year_label), 
           loadings_plot, width = 10, height = 8)
  }
  
  # Create summary table
  summary_table <- data.frame()
  for (system in names(pca_results_by_system)) {
    if (!is.null(pca_results_by_system[[system]])) {
      summary_table <- rbind(summary_table, data.frame(
        Year = year_label,
        System = system,
        N_Complete = pca_results_by_system[[system]]$n_complete,
        Correlation = round(abs(pca_results_by_system[[system]]$correlation), 3),
        Variance_PC1 = round(pca_results_by_system[[system]]$variance_explained[1] * 100, 1)
      ))
    }
  }
  
  return(list(
    results_by_system = pca_results_by_system,
    summary_table = summary_table
  ))
}

# Function to analyze asset ownership patterns by scoring system
analyze_assets_by_system <- function(data, year_label, scoring_analysis) {
  cat(sprintf("\n=== Asset Analysis by Scoring System for Year %s ===\n", year_label))
  
  data <- scoring_analysis$data
  existing_assets <- asset_vars[asset_vars %in% names(data)]
  
  # Calculate ownership rates by system
  ownership_by_system <- data.frame()
  
  for (system in unique(data$scoring_system_binary)) {
    system_data <- data %>% filter(scoring_system_binary == system)
    
    for (asset in existing_assets) {
      ownership_rate <- mean(system_data[[asset]] == 1, na.rm = TRUE) * 100
      ownership_by_system <- rbind(ownership_by_system, data.frame(
        Asset = asset,
        System = system,
        Ownership_Rate = ownership_rate
      ))
    }
  }
  
  # Create visualization
  ownership_plot <- ggplot(ownership_by_system, 
                           aes(x = reorder(Asset, Ownership_Rate), y = Ownership_Rate, 
                               fill = System)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(title = sprintf("Asset Ownership Rates by Scoring System - %s", year_label),
         x = "Asset",
         y = "Ownership Rate (%)") +
    theme_minimal()
  
  # Calculate differences (only if both systems exist)
  if (length(unique(ownership_by_system$System)) > 1) {
    ownership_diff <- ownership_by_system %>%
      pivot_wider(names_from = System, values_from = Ownership_Rate) %>%
      mutate(Difference = `System 1` - `System 2`) %>%
      arrange(desc(abs(Difference)))
  } else {
    # If only one system, create a simple summary
    ownership_diff <- ownership_by_system %>%
      pivot_wider(names_from = System, values_from = Ownership_Rate) %>%
      mutate(Difference = 0)  # No difference if only one system
  }
  
  ggsave(sprintf("descriptive_stats/asset_ownership_by_system_%s.pdf", year_label), 
         ownership_plot, width = 10, height = 8)
  
  return(list(
    ownership_by_system = ownership_by_system,
    ownership_diff = ownership_diff,
    plot = ownership_plot
  ))
}

# Main analysis pipeline
run_comprehensive_analysis <- function() {
  all_results <- list()
  
  for (year_name in c("y2000", "y2001", "y2002", "y2004")) {
    year_label <- gsub("y", "", year_name)
    cat(sprintf("\n\n========== ANALYZING YEAR %s ==========\n", year_label))
    
    # Step 1: Analyze scoring systems
    scoring_analysis <- analyze_scoring_systems(prepared_data[[year_name]], year_label)
    
    # Step 2: Perform PCA by scoring system
    pca_analysis <- perform_pca_by_system(prepared_data[[year_name]], year_label, scoring_analysis)
    
    # Step 3: Analyze assets by scoring system
    asset_analysis <- analyze_assets_by_system(prepared_data[[year_name]], year_label, scoring_analysis)
    
    # Store all results
    all_results[[year_name]] <- list(
      scoring = scoring_analysis,
      pca = pca_analysis,
      assets = asset_analysis
    )
    
    # Save detailed tables
    write.csv(scoring_analysis$system_summary, 
              sprintf("tables/scoring_system_summary_%s.csv", year_label), 
              row.names = FALSE)
    
    write.csv(scoring_analysis$regional_systems, 
              sprintf("tables/regional_scoring_systems_%s.csv", year_label), 
              row.names = FALSE)
    
    if (!is.null(pca_analysis$summary_table) && nrow(pca_analysis$summary_table) > 0) {
      write.csv(pca_analysis$summary_table, 
                sprintf("tables/pca_summary_by_system_%s.csv", year_label), 
                row.names = FALSE)
    }
    
    write.csv(asset_analysis$ownership_diff, 
              sprintf("tables/asset_ownership_differences_%s.csv", year_label), 
              row.names = FALSE)
  }
  
  # Create cross-year comparison
  create_cross_year_comparison(all_results)
  
  return(all_results)
}

# Function to create cross-year comparison
create_cross_year_comparison <- function(all_results) {
  cat("\n\n=== Creating Cross-Year Comparisons ===\n")
  
  # 1. Evolution of scoring systems
  system_evolution <- data.frame()
  for (year_name in names(all_results)) {
    year <- gsub("y", "", year_name)
    summary <- all_results[[year_name]]$scoring$system_summary
    if (nrow(summary) > 0) {
      summary$Year <- year
      system_evolution <- rbind(system_evolution, summary)
    }
  }
  
  # Plot evolution
  evolution_plot <- system_evolution %>%
    ggplot(aes(x = Year, y = pct, fill = scoring_system_binary, group = scoring_system_binary)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = "Evolution of Scoring Systems (2000-2004)",
         subtitle = "Percentage of households in each scoring system",
         y = "Percentage of Households",
         fill = "Scoring System") +
    theme_minimal()
  
  ggsave("pca_plots/scoring_system_evolution.pdf", evolution_plot, width = 10, height = 6)
  
  # 2. PCA correlation evolution
  pca_evolution <- data.frame()
  for (year_name in names(all_results)) {
    if (!is.null(all_results[[year_name]]$pca$summary_table)) {
      pca_evolution <- rbind(pca_evolution, all_results[[year_name]]$pca$summary_table)
    }
  }
  
  if (nrow(pca_evolution) > 0) {
    correlation_evolution_plot <- pca_evolution %>%
      ggplot(aes(x = Year, y = Correlation, color = System, group = System)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      labs(title = "PCA-PUN_INC Correlation Evolution by Scoring System",
           subtitle = "Higher correlation indicates stronger asset-based scoring",
           y = "Absolute Correlation Coefficient") +
      theme_minimal()
    
    ggsave("pca_plots/correlation_evolution_by_system.pdf", 
           correlation_evolution_plot, width = 10, height = 6)
  }
  
  # 3. Create summary report
  create_summary_report(all_results, system_evolution, pca_evolution)
}

# Function to create a summary report
create_summary_report <- function(all_results, system_evolution, pca_evolution) {
  report <- "# ENCASEH PCA Analysis Summary Report\n\n"
  report <- paste0(report, "## Key Findings\n\n")
  
  # Finding 1: Dual scoring system
  if (any(system_evolution$scoring_system_binary == "System 2")) {
    report <- paste0(report, "### 1. Dual Scoring System Identified\n")
    report <- paste0(report, "- Years 2000 and 2001 show evidence of two distinct scoring systems\n")
    report <- paste0(report, sprintf("- System 2 (high scores >100) was used for %.1f%% of households in 2000\n",
                                     system_evolution$pct[system_evolution$Year == "2000" & 
                                                            system_evolution$scoring_system_binary == "System 2"]))
    report <- paste0(report, "- By 2002, the dual system was eliminated\n\n")
  }
  
  # Finding 2: Low correlations
  report <- paste0(report, "### 2. PUN_INC is Not Purely Asset-Based\n")
  if (nrow(pca_evolution) > 0) {
    avg_correlation <- mean(pca_evolution$Correlation)
    report <- paste0(report, sprintf("- Average PCA-PUN_INC correlation: %.3f\n", avg_correlation))
    report <- paste0(report, "- Low correlations suggest PUN_INC incorporates factors beyond asset ownership\n\n")
  }
  
  # Finding 3: Regional patterns
  report <- paste0(report, "### 3. Regional Patterns in Scoring Systems\n")
  for (year_name in names(all_results)) {
    year <- gsub("y", "", year_name)
    n_system2_regions <- length(all_results[[year_name]]$scoring$system2_regions)
    if (n_system2_regions > 0) {
      report <- paste0(report, sprintf("- Year %s: %d regions predominantly used System 2\n", 
                                       year, n_system2_regions))
    }
  }
  
  # Save report
  writeLines(report, "tables/analysis_summary_report.md")
  cat("\nSummary report saved to tables/analysis_summary_report.md\n")
}

# Run the comprehensive analysis
all_results <- run_comprehensive_analysis()

# Create a final summary table for easy reference
final_summary <- data.frame(
  Metric = c("Dual Scoring System Present", 
             "System 2 Usage (%)",
             "Avg PCA Correlation (System 1)",
             "Avg PCA Correlation (System 2)",
             "Key Finding"),
  Year_2000 = c("Yes", 
                sprintf("%.1f", ifelse(any(all_results$y2000$scoring$system_summary$scoring_system_binary == "System 2"),
                                       all_results$y2000$scoring$system_summary$pct[all_results$y2000$scoring$system_summary$scoring_system_binary == "System 2"],
                                       0)),
                sprintf("%.3f", ifelse(any(names(all_results$y2000$pca$results_by_system) == "System 1"),
                                       abs(all_results$y2000$pca$results_by_system[["System 1"]]$correlation),
                                       NA)),
                sprintf("%.3f", ifelse(any(names(all_results$y2000$pca$results_by_system) == "System 2"),
                                       abs(all_results$y2000$pca$results_by_system[["System 2"]]$correlation),
                                       NA)),
                "Dual system active"),
  Year_2001 = c("Yes", 
                sprintf("%.1f", ifelse(any(all_results$y2001$scoring$system_summary$scoring_system_binary == "System 2"),
                                       all_results$y2001$scoring$system_summary$pct[all_results$y2001$scoring$system_summary$scoring_system_binary == "System 2"],
                                       0)),
                sprintf("%.3f", ifelse(any(names(all_results$y2001$pca$results_by_system) == "System 1"),
                                       abs(all_results$y2001$pca$results_by_system[["System 1"]]$correlation),
                                       NA)),
                sprintf("%.3f", ifelse(any(names(all_results$y2001$pca$results_by_system) == "System 2"),
                                       abs(all_results$y2001$pca$results_by_system[["System 2"]]$correlation),
                                       NA)),
                "Dual system declining"),
  Year_2002 = c("No", 
                "0.0",
                sprintf("%.3f", ifelse(any(names(all_results$y2002$pca$results_by_system) == "System 1"),
                                       abs(all_results$y2002$pca$results_by_system[["System 1"]]$correlation),
                                       NA)),
                "N/A",
                "Unified system"),
  Year_2004 = c("No", 
                "0.0",
                sprintf("%.3f", ifelse(any(names(all_results$y2004$pca$results_by_system) == "System 1"),
                                       abs(all_results$y2004$pca$results_by_system[["System 1"]]$correlation),
                                       NA)),
                "N/A",
                "Stable system")
)

write.csv(final_summary, "tables/final_summary_table.csv", row.names = FALSE)

cat("\n\nAnalysis complete! Check the following directories for outputs:\n")
cat("- descriptive_stats/: Distribution plots and asset analysis\n")
cat("- pca_plots/: PCA visualizations and correlations\n")
cat("- tables/: Detailed data tables and summary report\n")

# Create correlation summary visualization
create_correlation_summary_plot <- function(all_results) {
  
  # Extract correlations from all results
  correlation_data <- data.frame()
  
  for (year_name in names(all_results)) {
    year <- gsub("y", "", year_name)
    pca_results <- all_results[[year_name]]$pca$results_by_system
    
    if (!is.null(pca_results)) {
      for (system in names(pca_results)) {
        if (!is.null(pca_results[[system]])) {
          correlation_data <- rbind(correlation_data, data.frame(
            Year = as.numeric(year),
            System = system,
            Correlation = abs(pca_results[[system]]$correlation),  # Use absolute value
            Label = sprintf("%.3f", abs(pca_results[[system]]$correlation))
          ))
        }
      }
    }
  }
  
  # Create the plot
  correlation_plot <- ggplot(correlation_data, aes(x = factor(Year), y = Correlation, fill = System)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = Label), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, 
              size = 3.5) +
    scale_fill_manual(values = c("System 1" = "#2E86AB", "System 2" = "#E63946")) +
    labs(title = "PCA-PUN_INC Correlation by Year and Scoring System",
         subtitle = "Absolute correlation values shown\nSystem 1: Normal scores (<100), System 2: High scores (>100)",
         x = "Year",
         y = "Absolute Correlation",
         fill = "Scoring System") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10, color = "gray50"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    ylim(0, 0.35) +  # Set y-axis limit to show all values clearly
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50", alpha = 0.5) +
    geom_hline(yintercept = 0.2, linetype = "dashed", color = "gray50", alpha = 0.5)
  
  # Save the plot
  ggsave("pca_plots/pca_correlation_summary_by_system.pdf", 
         correlation_plot, 
         width = 10, 
         height = 6)
  
  # Also create a simple table version
  correlation_table <- correlation_data %>%
    select(-Label) %>%
    pivot_wider(names_from = System, values_from = Correlation) %>%
    arrange(Year)
  
  # Add a note about what's happening
  cat("\n=== PCA-PUN_INC Correlation Summary ===\n")
  print(correlation_table)
  cat("\nNote: Values are absolute correlations. Low values (<0.1) indicate PUN_INC is not primarily asset-based.\n")
  
  # Save the table
  write.csv(correlation_table, "pca_plots/correlation_summary_table.csv", row.names = FALSE)
  
  return(list(plot = correlation_plot, table = correlation_table))
}

# Call this function at the end of your run_comprehensive_analysis() function
# Right before the return statement, add:
correlation_summary <- create_correlation_summary_plot(all_results)

