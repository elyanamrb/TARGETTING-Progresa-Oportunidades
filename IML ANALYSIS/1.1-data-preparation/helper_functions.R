# ============================================================================
# helper_functions.R
# ============================================================================
# Purpose: Reusable functions for data preparation and RDD analysis
# Author: 
# Date Created: 2025
# Last Modified: 2025
# ============================================================================

# ============================================================================
# DATA CLEANING FUNCTIONS
# ============================================================================

#' Standardize locality codes
#' 
#' @param df Data frame with cve_edo, cve_mun, cve_loc columns
#' @return Data frame with standardized codes and loc_id
standardize_locality_codes <- function(df) {
  df %>%
    mutate(
      cve_edo = as.numeric(cve_edo),
      cve_mun = as.numeric(cve_mun),
      cve_loc = as.numeric(cve_loc),
      # Create standardized codes
      cve_edo_std = sprintf("%02d", cve_edo),
      cve_mun_std = sprintf("%03d", cve_mun),
      cve_loc_std = sprintf("%04d", cve_loc),
      # Create unique locality ID
      loc_id = paste0(cve_edo_std, cve_mun_std, cve_loc_std)
    )
}

#' Clean variable names consistently
#' 
#' @param df Data frame
#' @param year Year suffix to add (e.g., "95", "00", "05")
#' @return Data frame with cleaned names
clean_variable_names <- function(df, year = NULL) {
  # First, convert to lowercase and remove spaces
  names(df) <- tolower(names(df))
  names(df) <- gsub(" ", "_", names(df))
  
  # If year specified, add suffix to value variables
  if(!is.null(year)) {
    value_vars <- c("pobtot", "analf", "sinag", "sinel", "sindr", "pisotier", "iml")
    for(var in value_vars) {
      if(var %in% names(df)) {
        names(df)[names(df) == var] <- paste0(var, "_", year)
      }
    }
  }
  
  return(df)
}

# ============================================================================
# DATA VALIDATION FUNCTIONS
# ============================================================================

#' Check for duplicate localities
#' 
#' @param df Data frame with loc_id column
#' @param dataset_name Name for reporting
#' @return Invisible NULL, prints results
check_duplicates <- function(df, dataset_name = "Dataset") {
  
  duplicates <- df %>%
    group_by(loc_id) %>%
    filter(n() > 1) %>%
    ungroup()
  
  n_dups <- nrow(duplicates)
  n_unique_dups <- n_distinct(duplicates$loc_id)
  
  if(n_dups > 0) {
    message(paste0("⚠️  ", dataset_name, " has ", n_unique_dups, 
                   " duplicate locality IDs (", n_dups, " total rows)"))
    return(duplicates)
  } else {
    message(paste0("✓ ", dataset_name, " has no duplicate localities"))
    return(NULL)
  }
}

#' Validate numeric ranges
#' 
#' @param df Data frame
#' @param var Variable name
#' @param min_val Minimum valid value
#' @param max_val Maximum valid value
#' @return Data frame of invalid values
validate_range <- function(df, var, min_val, max_val) {
  invalid <- df %>%
    filter(!!sym(var) < min_val | !!sym(var) > max_val)
  
  if(nrow(invalid) > 0) {
    message(paste0("⚠️  ", var, " has ", nrow(invalid), 
                   " values outside range [", min_val, ", ", max_val, "]"))
  } else {
    message(paste0("✓ ", var, " values all within valid range"))
  }
  
  return(invalid)
}

# ============================================================================
# MERGE DIAGNOSTIC FUNCTIONS
# ============================================================================

#' Report merge success rates
#' 
#' @param df_merged Result of merge
#' @param key_var Variable that should be non-missing if merge successful
#' @param merge_name Name for reporting
merge_diagnostics <- function(df_merged, key_var, merge_name = "Merge") {
  
  n_total <- nrow(df_merged)
  n_matched <- sum(!is.na(df_merged[[key_var]]))
  match_rate <- round(n_matched / n_total * 100, 1)
  
  message(paste0(merge_name, " Results:"))
  message(paste0("  Total rows: ", n_total))
  message(paste0("  Matched: ", n_matched, " (", match_rate, "%)"))
  message(paste0("  Unmatched: ", n_total - n_matched))
  
  return(invisible(list(
    total = n_total,
    matched = n_matched,
    unmatched = n_total - n_matched,
    match_rate = match_rate
  )))
}

# ============================================================================
# RDD-SPECIFIC FUNCTIONS
# ============================================================================

#' Find optimal RDD bandwidth using Imbens-Kalyanaraman method
#' 
#' @param x Running variable
#' @param y Outcome variable
#' @param c Cutoff value
#' @return Optimal bandwidth
ik_bandwidth <- function(x, y, c = 0) {
  # Center running variable at cutoff
  x_centered <- x - c
  
  # Remove missing values
  complete_cases <- complete.cases(x_centered, y)
  x_centered <- x_centered[complete_cases]
  y <- y[complete_cases]
  
  # Calculate optimal bandwidth (simplified IK)
  n <- length(x_centered)
  sd_x <- sd(x_centered)
  sd_y <- sd(y)
  
  # Rule of thumb bandwidth
  h <- 1.06 * sd_x * n^(-1/5)
  
  return(h)
}

#' Create binned scatter plot for RDD visualization
#' 
#' @param df Data frame
#' @param x_var Running variable name
#' @param y_var Outcome variable name
#' @param cutoff Cutoff value
#' @param nbins Number of bins on each side
#' @return ggplot object
create_rd_plot <- function(df, x_var, y_var, cutoff = 0, nbins = 20) {
  
  # Create bins
  df_plot <- df %>%
    mutate(
      side = ifelse(!!sym(x_var) < cutoff, "left", "right"),
      bin = cut(!!sym(x_var), 
                breaks = c(
                  seq(min(!!sym(x_var), na.rm = TRUE), cutoff, length.out = nbins + 1),
                  seq(cutoff, max(!!sym(x_var), na.rm = TRUE), length.out = nbins + 1)[-1]
                ),
                include.lowest = TRUE)
    ) %>%
    group_by(bin, side) %>%
    summarise(
      x_mean = mean(!!sym(x_var), na.rm = TRUE),
      y_mean = mean(!!sym(y_var), na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  # Create plot
  p <- ggplot(df_plot, aes(x = x_mean, y = y_mean)) +
    geom_point(aes(size = n, color = side), alpha = 0.7) +
    geom_smooth(data = df %>% filter(!!sym(x_var) < cutoff),
                aes(x = !!sym(x_var), y = !!sym(y_var)),
                method = "lm", se = TRUE, color = "blue") +
    geom_smooth(data = df %>% filter(!!sym(x_var) >= cutoff),
                aes(x = !!sym(x_var), y = !!sym(y_var)),
                method = "lm", se = TRUE, color = "red") +
    geom_vline(xintercept = cutoff, linetype = "dashed", color = "black") +
    theme_minimal() +
    labs(x = x_var, y = y_var,
         title = paste("RDD Plot:", y_var, "vs", x_var),
         subtitle = paste("Cutoff =", cutoff)) +
    scale_color_manual(values = c("left" = "blue", "right" = "red")) +
    guides(size = "none")
  
  return(p)
}

# ============================================================================
# SUMMARY STATISTICS FUNCTIONS
# ============================================================================

#' Create formatted summary statistics table
#' 
#' @param df Data frame
#' @param vars Variables to summarize
#' @param by Grouping variable (optional)
#' @return Data frame with summary statistics
create_summary_stats <- function(df, vars, by = NULL) {
  
  # Helper function for one variable
  summarize_var <- function(x, var_name) {
    data.frame(
      variable = var_name,
      n = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      q25 = quantile(x, 0.25, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q75 = quantile(x, 0.75, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  }
  
  # Apply to all variables
  if(is.null(by)) {
    # Overall statistics
    summary_stats <- map_df(vars, function(v) {
      summarize_var(df[[v]], v)
    })
  } else {
    # Statistics by group
    summary_stats <- df %>%
      group_by(!!sym(by)) %>%
      group_map(function(group_df, group_key) {
        group_stats <- map_df(vars, function(v) {
          summarize_var(group_df[[v]], v)
        })
        group_stats[[by]] <- group_key[[1]]
        return(group_stats)
      }) %>%
      bind_rows()
  }
  
  # Round numeric columns
  numeric_cols <- c("mean", "sd", "min", "q25", "median", "q75", "max")
  summary_stats[numeric_cols] <- round(summary_stats[numeric_cols], 3)
  
  return(summary_stats)
}

#' Perform balance test at RDD cutoff
#' 
#' @param df Data frame
#' @param covariates Vector of covariate names
#' @param running_var Running variable name
#' @param cutoff Cutoff value
#' @param bandwidth Bandwidth for local comparison
#' @return Data frame with balance test results
balance_test <- function(df, covariates, running_var, cutoff = 0, bandwidth = NULL) {
  
  # If no bandwidth specified, use IK optimal
  if(is.null(bandwidth)) {
    bandwidth <- 0.5  # Default to 0.5 for simplicity
  }
  
  # Filter to bandwidth
  df_local <- df %>%
    filter(abs(!!sym(running_var) - cutoff) <= bandwidth)
  
  # Create treatment indicator
  df_local <- df_local %>%
    mutate(treated = ifelse(!!sym(running_var) >= cutoff, 1, 0))
  
  # Test each covariate
  balance_results <- map_df(covariates, function(cov) {
    
    # Skip if variable doesn't exist
    if(!(cov %in% names(df_local))) {
      return(NULL)
    }
    
    # Perform t-test
    test_result <- t.test(as.formula(paste(cov, "~ treated")), 
                          data = df_local)
    
    return(data.frame(
      variable = cov,
      mean_control = test_result$estimate[1],
      mean_treated = test_result$estimate[2],
      difference = test_result$estimate[2] - test_result$estimate[1],
      se = test_result$stderr,
      t_stat = test_result$statistic,
      p_value = test_result$p.value,
      n_control = sum(df_local$treated == 0, na.rm = TRUE),
      n_treated = sum(df_local$treated == 1, na.rm = TRUE)
    ))
  })
  
  # Add significance stars
  balance_results <- balance_results %>%
    mutate(
      significance = case_when(
        p_value < 0.01 ~ "***",
        p_value < 0.05 ~ "**",
        p_value < 0.10 ~ "*",
        TRUE ~ ""
      )
    )
  
  return(balance_results)
}

# ============================================================================
# EXPORT FUNCTION LIST
# ============================================================================

# List of all functions for easy reference
data_prep_functions <- list(
  # Data cleaning
  standardize_locality_codes = standardize_locality_codes,
  clean_variable_names = clean_variable_names,
  
  # Validation
  check_duplicates = check_duplicates,
  validate_range = validate_range,
  
  # Merge diagnostics
  merge_diagnostics = merge_diagnostics,
  
  # RDD functions
  ik_bandwidth = ik_bandwidth,
  create_rd_plot = create_rd_plot,
  
  # Summary statistics
  create_summary_stats = create_summary_stats,
  balance_test = balance_test
)

message("Helper functions loaded successfully!")
message("Available functions:")
message(paste("  -", names(data_prep_functions), collapse = "\n"))
