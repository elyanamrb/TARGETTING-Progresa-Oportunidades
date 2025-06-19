# 03_ventanilla_rdd_analysis.R
# RDD Analysis using Ventanilla (self-selection) as treatment variable
# Focus on 2004 data only

# Load required libraries
library(tidyverse)
library(rdrobust)
library(stargazer)
library(gridExtra)
library(scales)

# Set working directory
setwd("/Users/elyanaramos/Desktop/TA:RA/RA/V4")

# Load prepared data
prepared_data <- readRDS("prepared_encaseh_data.rds")

# Extract 2004 data
data_2004 <- prepared_data$y2004

# ========================================
# PART 1: DATA PREPARATION
# ========================================

# Create ventanilla dummy variable (1 = Ventanilla, 0 = Other methods)
data_2004 <- data_2004 %>%
  filter(!is.na(targeting_method)) %>%  # Remove NA targeting methods
  mutate(
    ventanilla = ifelse(targeting_method == "Ventanilla", 1, 0),
    # Clean PUN_INC to remove extreme outliers for cleaner analysis
    pun_inc_clean = ifelse(PUN_INC > 100, NA, PUN_INC)
  )

# Print summary statistics
cat("\n========== 2004 Data Summary ==========\n")
cat("Total observations with targeting method:", nrow(data_2004), "\n")
cat("Ventanilla observations:", sum(data_2004$ventanilla), "\n")
cat("Other methods observations:", sum(1 - data_2004$ventanilla), "\n")
cat("Proportion using Ventanilla:", mean(data_2004$ventanilla), "\n\n")

# ========================================
# PART 2: NATIONAL RDD ANALYSIS
# ========================================

# Function to find optimal cutoff for RDD - OPTIMIZED VERSION
find_optimal_cutoff <- function(data, outcome_var = "ventanilla", running_var = "pun_inc_clean",
                                cutoff_range = seq(-5, 5, by = 0.1), sample_size = 50000) {
  
  # Remove NA values
  analysis_data <- data %>%
    filter(!is.na(!!sym(outcome_var)) & !is.na(!!sym(running_var)))
  
  # Skip if insufficient data
  if (nrow(analysis_data) < 50) return(NULL)
  
  # OPTIMIZATION: Take a random sample if data is too large
  if (nrow(analysis_data) > sample_size) {
    cat("  Data has", nrow(analysis_data), "observations. Sampling", sample_size, "for cutoff search...\n")
    set.seed(42)  # For reproducibility
    sample_data <- analysis_data %>%
      sample_n(sample_size)
  } else {
    sample_data <- analysis_data
  }
  
  # First, get a rough idea of where the cutoff might be by examining the data
  # Look at the distribution of the running variable by treatment
  treated_mean <- mean(sample_data[[running_var]][sample_data[[outcome_var]] == 1], na.rm = TRUE)
  control_mean <- mean(sample_data[[running_var]][sample_data[[outcome_var]] == 0], na.rm = TRUE)
  
  # Create a more focused cutoff range around the means
  rough_center <- (treated_mean + control_mean) / 2
  focused_range <- seq(rough_center - 2, rough_center + 2, by = 0.2)
  
  cat("  Searching for optimal cutoff in range:", min(focused_range), "to", max(focused_range), "\n")
  
  # Initialize results storage
  results <- data.frame()
  
  # Loop over possible cutoffs
  max_effect <- -Inf
  optimal_cutoff <- NA
  optimal_result <- NULL
  
  # Progress counter
  n_cutoffs <- length(focused_range)
  counter <- 0
  
  for (cutoff in focused_range) {
    counter <- counter + 1
    if (counter %% 5 == 0) {
      cat("  Progress:", round(counter/n_cutoffs * 100), "%\n")
    }
    
    tryCatch({
      # Perform RD analysis with smaller bandwidth for speed
      rd_result <- rdrobust(
        y = sample_data[[outcome_var]],
        x = sample_data[[running_var]],
        c = cutoff,
        kernel = "triangular",
        bwselect = "mserd",
        masspoints = "off",  # Turn off mass points check for speed
        all = FALSE  # Don't compute all statistics
      )
      
      # Store results
      effect_size <- abs(rd_result$coef[1])
      results <- rbind(results, data.frame(
        cutoff = cutoff,
        effect = rd_result$coef[1],
        se = rd_result$se[1],
        pvalue = rd_result$pv[1],
        n_left = rd_result$N_h[1],
        n_right = rd_result$N_h[2]
      ))
      
      # Update optimal cutoff if this one has larger effect and is significant
      if (!is.na(effect_size) && effect_size > max_effect && rd_result$pv[1] < 0.1) {
        max_effect <- effect_size
        optimal_cutoff <- cutoff
      }
      
    }, error = function(e) NULL)
  }
  
  # If no significant results found, use cutoff with maximum effect
  if (is.na(optimal_cutoff) && nrow(results) > 0) {
    max_idx <- which.max(abs(results$effect))
    optimal_cutoff <- results$cutoff[max_idx]
  }
  
  cat("  Optimal cutoff found:", optimal_cutoff, "\n")
  
  # Now run the final RD with the full dataset at the optimal cutoff
  cat("  Running final RD with full dataset...\n")
  optimal_result <- rdrobust(
    y = analysis_data[[outcome_var]],
    x = analysis_data[[running_var]],
    c = optimal_cutoff,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  return(list(
    optimal_cutoff = optimal_cutoff,
    rd_result = optimal_result,
    all_results = results,
    data = analysis_data,
    sample_size_used = nrow(sample_data)
  ))
}

# Alternative: Quick single cutoff analysis
quick_rdd_analysis <- function(data, outcome_var = "ventanilla", running_var = "pun_inc_clean",
                               cutoff = NULL) {
  
  # Remove NA values
  analysis_data <- data %>%
    filter(!is.na(!!sym(outcome_var)) & !is.na(!!sym(running_var)))
  
  # If no cutoff provided, use median of running variable
  if (is.null(cutoff)) {
    cutoff <- median(analysis_data[[running_var]], na.rm = TRUE)
    cat("  Using median as cutoff:", cutoff, "\n")
  }
  
  # Run RDD
  rd_result <- rdrobust(
    y = analysis_data[[outcome_var]],
    x = analysis_data[[running_var]],
    c = cutoff,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  return(list(
    cutoff = cutoff,
    rd_result = rd_result,
    data = analysis_data
  ))
}

# Run national RDD analysis
cat("Finding optimal cutoff for national analysis...\n")

# Option 1: Use the optimized search (recommended)
national_results <- find_optimal_cutoff(data_2004, sample_size = 100000)

# Option 2: If still too slow, use quick analysis with pre-specified cutoff
# You can use the median or a theoretically motivated cutoff
# national_results <- quick_rdd_analysis(data_2004, cutoff = 0)

if (!is.null(national_results)) {
  cat("\nNational RDD Results:\n")
  cat("Optimal cutoff:", round(national_results$optimal_cutoff, 3), "\n")
  cat("RD Estimate:", round(national_results$rd_result$coef[1], 3), "\n")
  cat("Standard Error:", round(national_results$rd_result$se[1], 3), "\n")
  cat("P-value:", round(national_results$rd_result$pv[1], 3), "\n")
  cat("N observations:", nrow(national_results$data), "\n")
  if (exists("sample_size_used", where = national_results)) {
    cat("Sample size used for search:", national_results$sample_size_used, "\n")
  }
  cat("\n")
}

# ========================================
# PART 3: VISUALIZATION FUNCTIONS
# ========================================

create_rdd_plot <- function(data, cutoff, outcome_var = "ventanilla", 
                            running_var = "pun_inc_clean", title_prefix = "National") {
  
  # Create binned data for visualization
  # Use more bins for smoother visualization
  n_bins <- 40
  bin_width <- (max(data[[running_var]], na.rm = TRUE) - min(data[[running_var]], na.rm = TRUE)) / n_bins
  
  plot_data <- data %>%
    mutate(bin = cut(!!sym(running_var), 
                     breaks = seq(min(!!sym(running_var), na.rm = TRUE), 
                                  max(!!sym(running_var), na.rm = TRUE), 
                                  length.out = n_bins + 1),
                     include.lowest = TRUE))
  
  bin_means <- plot_data %>%
    group_by(bin) %>%
    summarize(
      mean_x = mean(!!sym(running_var), na.rm = TRUE),
      mean_y = mean(!!sym(outcome_var), na.rm = TRUE),
      se_y = sd(!!sym(outcome_var), na.rm = TRUE) / sqrt(n()),
      n = n(),
      .groups = "drop"
    ) %>%
    filter(n >= 10)  # Require at least 10 observations per bin
  
  # Determine significance for color coding
  rd_test <- rdrobust(
    y = data[[outcome_var]],
    x = data[[running_var]],
    c = cutoff,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  cutoff_color <- ifelse(rd_test$pv[1] < 0.05, "red", "gray50")
  significance_label <- ifelse(rd_test$pv[1] < 0.05, "Significant (p < 0.05)", "Not significant")
  
  # Create RDD plot with proper y-axis limits
  p <- ggplot(bin_means, aes(x = mean_x, y = mean_y)) +
    geom_point(color = "steelblue", size = 3, alpha = 0.7) +
    geom_errorbar(aes(ymin = mean_y - se_y, ymax = mean_y + se_y), 
                  width = bin_width/2, alpha = 0.3) +
    geom_vline(xintercept = cutoff, linetype = "dashed", color = cutoff_color, size = 1) +
    geom_smooth(data = subset(bin_means, mean_x < cutoff),
                method = "lm", se = TRUE, color = "cornflowerblue") +
    geom_smooth(data = subset(bin_means, mean_x >= cutoff),
                method = "lm", se = TRUE, color = "coral") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +  # Force y-axis to be 0-1
    labs(
      title = paste(title_prefix, "- Ventanilla RDD (2004)"),
      subtitle = paste("Cutoff:", round(cutoff, 3), 
                       "| Effect:", round(rd_test$coef[1], 3),
                       "| SE:", round(rd_test$se[1], 3),
                       "| ", significance_label),
      x = "Poverty Score (PUN_INC)",
      y = "Probability of Using Ventanilla"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11)
    )
  
  return(p)
}

# Create national RDD plot
if (!is.null(national_results)) {
  national_plot <- create_rdd_plot(
    national_results$data, 
    national_results$optimal_cutoff,
    title_prefix = "National"
  )
  
  # Save plot
  ggsave("output/national_ventanilla_rdd_2004.pdf", national_plot, 
         width = 10, height = 6)
}

# ========================================
# PART 4: STATE-LEVEL ANALYSIS
# ========================================

# Get state names for better labeling
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

# Create output directory for state results
dir.create("output/state_rdd_results", showWarnings = FALSE, recursive = TRUE)

# Initialize state results storage
state_results_table <- data.frame()

# Get unique states
unique_states <- sort(unique(data_2004$ENT_FED))

# Run analysis for each state
cat("\nRunning state-level RDD analysis...\n")
for (state_id in unique_states) {
  
  # Get state name
  state_name <- state_names[as.character(state_id)]
  if (is.na(state_name)) state_name <- paste("State", state_id)
  
  cat("Analyzing", state_name, "(", state_id, ")...\n")
  
  # Filter data for this state
  state_data <- data_2004 %>%
    filter(ENT_FED == state_id)
  
  # Check if state has both ventanilla and other methods
  n_ventanilla <- sum(state_data$ventanilla)
  n_other <- sum(1 - state_data$ventanilla)
  
  if (n_ventanilla < 10 || n_other < 10) {
    cat("  Insufficient variation in targeting methods for state", state_id, "\n")
    next
  }
  
  # Run RDD for this state
  state_results <- find_optimal_cutoff(state_data)
  
  if (is.null(state_results)) {
    cat("  No valid results for state", state_id, "\n")
    next
  }
  
  # Create and save state plot
  state_plot <- create_rdd_plot(
    state_results$data,
    state_results$optimal_cutoff,
    title_prefix = state_name
  )
  
  state_str <- sprintf("%02d", state_id)
  ggsave(paste0("output/state_rdd_results/state_", state_str, "_", 
                gsub(" ", "_", state_name), "_ventanilla_rdd.pdf"), 
         state_plot, width = 10, height = 6)
  
  # Add to results table
  state_results_table <- rbind(state_results_table,
                               data.frame(
                                 state_id = state_id,
                                 state_name = state_name,
                                 optimal_cutoff = state_results$optimal_cutoff,
                                 n_obs = nrow(state_results$data),
                                 n_ventanilla = sum(state_results$data$ventanilla),
                                 n_other = sum(1 - state_results$data$ventanilla),
                                 rd_estimate = state_results$rd_result$coef[1],
                                 std_error = state_results$rd_result$se[1],
                                 p_value = state_results$rd_result$pv[1],
                                 significant = state_results$rd_result$pv[1] < 0.05
                               ))
}

# ========================================
# PART 5: CREATE SUMMARY TABLES AND PLOTS
# ========================================

# Format state results table
if (nrow(state_results_table) > 0) {
  formatted_state_table <- state_results_table %>%
    mutate(
      rd_estimate = round(rd_estimate, 3),
      std_error = round(std_error, 3),
      p_value = round(p_value, 3),
      significance = ifelse(p_value < 0.01, "***", 
                            ifelse(p_value < 0.05, "**",
                                   ifelse(p_value < 0.1, "*", ""))),
      result = paste0(rd_estimate, significance),
      optimal_cutoff = round(optimal_cutoff, 2)
    ) %>%
    select(state_id, state_name, optimal_cutoff, n_obs, n_ventanilla, 
           result, std_error, p_value)
  
  # Save as CSV
  write.csv(formatted_state_table, "output/state_ventanilla_rdd_results.csv", 
            row.names = FALSE)
  
  # Generate LaTeX table
  stargazer(
    formatted_state_table %>% select(-state_id),
    type = "latex",
    summary = FALSE,
    rownames = FALSE,
    title = "State-Level Ventanilla RDD Results (2004)",
    digits = 3,
    align = TRUE,
    out = "output/state_ventanilla_rdd_results.tex"
  )
}

# ========================================
# PART 6: COMPARATIVE VISUALIZATION
# ========================================

# Create coefficient plot comparing significant states
if (nrow(state_results_table) > 0) {
  
  # Filter for states with sufficient observations
  plot_data <- state_results_table %>%
    filter(n_obs >= 100) %>%
    arrange(rd_estimate)
  
  coef_plot <- ggplot(plot_data, aes(x = reorder(state_name, rd_estimate), 
                                     y = rd_estimate)) +
    geom_point(aes(color = significant), size = 3) +
    geom_errorbar(aes(ymin = rd_estimate - 1.96*std_error,
                      ymax = rd_estimate + 1.96*std_error),
                  width = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "coral"),
                       labels = c("Not Significant", "Significant (p < 0.05)")) +
    coord_flip() +
    labs(
      title = "State-Level RDD Estimates: Effect of Poverty Score on Ventanilla Use",
      subtitle = "2004 Data - States with n ≥ 100",
      x = "",
      y = "RD Estimate (with 95% CI)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 14, face = "bold")
    )
  
  ggsave("output/state_coefficients_comparison.pdf", coef_plot, 
         width = 10, height = 8)
}

# ========================================
# PART 7: DISTRIBUTION ANALYSIS
# ========================================

# Create distribution plots comparing PUN_INC by targeting method
dist_plot_national <- ggplot(data_2004 %>% filter(!is.na(pun_inc_clean)), 
                             aes(x = pun_inc_clean, fill = factor(ventanilla))) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("0" = "steelblue", "1" = "coral"),
                    labels = c("Other Methods", "Ventanilla")) +
  labs(
    title = "Distribution of Poverty Scores by Targeting Method",
    subtitle = "National Level - 2004",
    x = "Poverty Score (PUN_INC)",
    y = "Density",
    fill = "Targeting Method"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("output/pun_inc_distribution_by_method_2004.pdf", dist_plot_national,
       width = 10, height = 6)

# ========================================
# PART 8: SUMMARY STATISTICS
# ========================================

# Create summary statistics table
summary_stats <- data_2004 %>%
  filter(!is.na(pun_inc_clean)) %>%
  group_by(ventanilla) %>%
  summarise(
    n = n(),
    mean_score = mean(pun_inc_clean, na.rm = TRUE),
    sd_score = sd(pun_inc_clean, na.rm = TRUE),
    median_score = median(pun_inc_clean, na.rm = TRUE),
    q25 = quantile(pun_inc_clean, 0.25, na.rm = TRUE),
    q75 = quantile(pun_inc_clean, 0.75, na.rm = TRUE),
    acceptance_rate = mean(accepted, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    targeting_method = ifelse(ventanilla == 1, "Ventanilla", "Other Methods")
  ) %>%
  select(targeting_method, everything(), -ventanilla)

# Print summary
cat("\n========== Summary Statistics by Targeting Method ==========\n")
print(summary_stats)

# Save summary statistics
write.csv(summary_stats, "output/summary_stats_by_method_2004.csv", row.names = FALSE)

# ========================================
# PART 9: ROBUSTNESS CHECKS
# ========================================

# Alternative RDD with RESULTADO as outcome
cat("\n\nRunning robustness check with RESULTADO as outcome...\n")

# Create alternative outcome analysis function with error handling
analyze_resultado_rdd <- function(data, cutoff, min_obs_around_cutoff = 50) {
  # Filter for valid data
  analysis_data <- data %>%
    filter(!is.na(accepted) & !is.na(pun_inc_clean))
  
  # Check if there's enough data around the cutoff
  left_obs <- sum(analysis_data$pun_inc_clean < cutoff)
  right_obs <- sum(analysis_data$pun_inc_clean >= cutoff)
  
  if (left_obs < min_obs_around_cutoff || right_obs < min_obs_around_cutoff) {
    cat("  Insufficient observations around cutoff (left:", left_obs, ", right:", right_obs, ")\n")
    return(NULL)
  }
  
  # Run RDD with error handling
  rd_result <- tryCatch({
    rdrobust(
      y = analysis_data$accepted,
      x = analysis_data$pun_inc_clean,
      c = cutoff,
      kernel = "triangular",
      bwselect = "mserd",
      masspoints = "off"  # Turn off mass points check to avoid errors
    )
  }, error = function(e) {
    cat("  Error in RDD estimation:", e$message, "\n")
    return(NULL)
  })
  
  return(rd_result)
}

# Run for ventanilla vs other methods
if (!is.null(national_results)) {
  # Ventanilla subset
  ventanilla_data <- data_2004 %>% filter(ventanilla == 1)
  other_data <- data_2004 %>% filter(ventanilla == 0)
  
  cat("\n=== Robustness Check: RESULTADO as Outcome ===\n")
  cat("Using cutoff:", round(national_results$optimal_cutoff, 3), "\n")
  
  # Check data availability
  cat("\nData summary:\n")
  cat("  Ventanilla observations:", nrow(ventanilla_data), "\n")
  cat("  Other methods observations:", nrow(other_data), "\n")
  
  # Run RDD for each group using the national cutoff
  cat("\nVentanilla Method:\n")
  ventanilla_resultado <- analyze_resultado_rdd(ventanilla_data, national_results$optimal_cutoff)
  
  if (!is.null(ventanilla_resultado)) {
    cat("  RD Estimate:", round(ventanilla_resultado$coef[1], 3), "\n")
    cat("  P-value:", round(ventanilla_resultado$pv[1], 3), "\n")
  } else {
    cat("  Could not estimate RDD for Ventanilla subset\n")
  }
  
  cat("\nOther Methods:\n")
  other_resultado <- analyze_resultado_rdd(other_data, national_results$optimal_cutoff)
  
  if (!is.null(other_resultado)) {
    cat("  RD Estimate:", round(other_resultado$coef[1], 3), "\n")
    cat("  P-value:", round(other_resultado$pv[1], 3), "\n")
  } else {
    cat("  Could not estimate RDD for Other Methods subset\n")
  }
  
  # Alternative: Compare acceptance rates by targeting method at the cutoff
  cat("\n=== Alternative Analysis: Acceptance Rates by Method ===\n")
  
  # Define bandwidth around cutoff for comparison
  bw <- 0.5  # You can adjust this
  
  comparison_data <- data_2004 %>%
    filter(!is.na(pun_inc_clean) & !is.na(accepted)) %>%
    filter(abs(pun_inc_clean - national_results$optimal_cutoff) <= bw) %>%
    mutate(
      side = ifelse(pun_inc_clean < national_results$optimal_cutoff, "Left", "Right"),
      method = ifelse(ventanilla == 1, "Ventanilla", "Other")
    )
  
  acceptance_summary <- comparison_data %>%
    group_by(method, side) %>%
    summarise(
      n = n(),
      acceptance_rate = mean(accepted, na.rm = TRUE),
      se = sd(accepted, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  print(acceptance_summary)
  
  # Create a simple comparison plot
  if (nrow(acceptance_summary) > 0) {
    comparison_plot <- ggplot(acceptance_summary, aes(x = side, y = acceptance_rate, 
                                                      color = method, group = method)) +
      geom_point(size = 4) +
      geom_line(size = 1.2) +
      geom_errorbar(aes(ymin = acceptance_rate - 1.96*se, 
                        ymax = acceptance_rate + 1.96*se), width = 0.1) +
      scale_color_manual(values = c("Ventanilla" = "coral", "Other" = "steelblue")) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(
        title = "Acceptance Rates by Targeting Method Around Cutoff",
        subtitle = paste("Within", bw, "of cutoff =", round(national_results$optimal_cutoff, 3)),
        x = "Side of Cutoff",
        y = "Acceptance Rate",
        color = "Targeting Method"
      ) +
      theme_minimal()
    
    ggsave("output/acceptance_rates_comparison.pdf", comparison_plot, 
           width = 8, height = 6)
  }
}

cat("\n\nAnalysis complete! All results saved in the 'output' folder.\n")
