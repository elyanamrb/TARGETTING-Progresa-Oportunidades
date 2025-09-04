# ============================================================================
# 02_robustness_national_rdd.R
# ============================================================================
# Purpose: Robustness checks for national RDD analysis
# Author: Elyana Ramos
# Date Created: 2025
# Last Modified: 2025
#
# Robustness tests:
#   1. Donut hole RDD (exclude observations near cutoff)
#   2. Alternative bandwidths
#   3. Placebo tests with different years
#   4. Alternative kernels
#   5. Higher-order polynomials
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(rdrobust)
library(gridExtra)

# ============================================================================
# SETUP
# ============================================================================

# Load the data and regional cutoffs
merged_entry_iml_rg <- read.csv("~/Desktop/TARGETTING-Progresa-Oportunidades/IML ANALYSIS/1.1-data-preparation/data/processed/merged_entry_iml_rg.csv")

# Function to get regional cutoffs (from main analysis)
get_regional_cutoffs <- function(data, post_year = 2002) {
  regions <- sort(unique(data$region41))
  regional_cutoffs <- data.frame()
  
  for (region_num in regions) {
    region_data <- data %>% 
      filter(region41 == region_num) %>%
      filter(!is.na(iml_95) & !is.na(entry_year_30))
    
    if (nrow(region_data) < 30) next
    
    region_data$post_year <- ifelse(region_data$entry_year_30 >= post_year, 1, 0)
    
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
  return(regional_cutoffs)
}

# Get cutoffs and normalize data
regional_cutoffs <- get_regional_cutoffs(merged_entry_iml_rg, 2002)
normalized_data <- merged_entry_iml_rg %>%
  mutate(post_2002 = ifelse(entry_year_30 >= 2002, 1, 0)) %>%
  inner_join(regional_cutoffs, by = c("region41" = "region")) %>%
  mutate(normalized_iml = iml_95 - cutoff) %>%
  filter(!is.na(normalized_iml),
         !is.na(entry_year_30),
         entry_year_30 >= 1995 & entry_year_30 <= 2010)

cat("========================================\n")
cat("NATIONAL RDD ROBUSTNESS CHECKS\n")
cat("========================================\n\n")

# Store main result for comparison
main_result <- rdrobust(
  y = normalized_data$post_2002,
  x = normalized_data$normalized_iml,
  c = 0,
  kernel = "triangular",
  bwselect = "mserd"
)

cat("MAIN RESULT (for comparison):\n")
cat("Estimate:", round(main_result$coef[1], 3), "\n")
cat("SE:", round(main_result$se[1], 3), "\n")
cat("P-value:", round(main_result$pv[1], 4), "\n\n")

# ============================================================================
# TEST 1: DONUT HOLE RDD
# ============================================================================

cat("------------------------------------\n")
cat("TEST 1: DONUT HOLE RDD\n")
cat("------------------------------------\n")
cat("Excluding observations near cutoff\n\n")

donut_sizes <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25)
donut_results <- data.frame()

for (donut in donut_sizes) {
  # Filter out observations within donut hole
  if (donut > 0) {
    donut_data <- normalized_data %>%
      filter(abs(normalized_iml) >= donut)
  } else {
    donut_data <- normalized_data
  }
  
  # Run RDD
  if (nrow(donut_data) >= 100) {
    rd_result <- tryCatch({
      rdrobust(
        y = donut_data$post_2002,
        x = donut_data$normalized_iml,
        c = 0,
        kernel = "triangular",
        bwselect = "mserd"
      )
    }, error = function(e) NULL)
    
    if (!is.null(rd_result)) {
      donut_results <- rbind(donut_results, data.frame(
        donut_size = donut,
        estimate = rd_result$coef[1],
        se = rd_result$se[1],
        pvalue = rd_result$pv[1],
        n_obs = nrow(donut_data),
        n_effective = rd_result$N_h[1] + rd_result$N_h[2]
      ))
      
      cat(sprintf("Donut = %.2f: Est = %.3f (SE = %.3f), p = %.4f, N = %d\n",
                  donut, rd_result$coef[1], rd_result$se[1], 
                  rd_result$pv[1], nrow(donut_data)))
    }
  }
}

# Plot donut results
if (nrow(donut_results) > 0) {
  p_donut <- ggplot(donut_results, aes(x = donut_size, y = estimate)) +
    geom_line(color = "blue") +
    geom_point(size = 3) +
    geom_ribbon(aes(ymin = estimate - 1.96*se, ymax = estimate + 1.96*se), 
                alpha = 0.2, fill = "blue") +
    geom_hline(yintercept = main_result$coef[1], linetype = "dashed", color = "red") +
    labs(title = "Donut Hole RDD: Sensitivity to Excluding Near-Cutoff Observations",
         x = "Donut Hole Size",
         y = "RD Estimate",
         caption = "Red dashed line = main estimate") +
    theme_minimal()
  
  ggsave("output/robustness_donut_national.pdf", p_donut, width = 10, height = 6)
}

# ============================================================================
# TEST 2: ALTERNATIVE BANDWIDTHS
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 2: ALTERNATIVE BANDWIDTHS\n")
cat("------------------------------------\n\n")

# Test different bandwidth selectors
bw_methods <- c("mserd", "msetwo", "msesum", "msecomb1", "msecomb2", "cerrd", "certwo")
bw_results <- data.frame()

for (method in bw_methods) {
  rd_result <- tryCatch({
    rdrobust(
      y = normalized_data$post_2002,
      x = normalized_data$normalized_iml,
      c = 0,
      kernel = "triangular",
      bwselect = method
    )
  }, error = function(e) NULL)
  
  if (!is.null(rd_result)) {
    bw_results <- rbind(bw_results, data.frame(
      method = method,
      bandwidth = rd_result$bws[1],
      estimate = rd_result$coef[1],
      se = rd_result$se[1],
      pvalue = rd_result$pv[1]
    ))
    
    cat(sprintf("%-10s: BW = %.3f, Est = %.3f (SE = %.3f), p = %.4f\n",
                method, rd_result$bws[1], rd_result$coef[1], 
                rd_result$se[1], rd_result$pv[1]))
  }
}

# ============================================================================
# TEST 3: PLACEBO TESTS (DIFFERENT YEARS)
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 3: PLACEBO TESTS\n")
cat("------------------------------------\n")
cat("Testing discontinuities at different years\n\n")

years <- 2000:2006
placebo_results <- data.frame()

for (year in years) {
  # Create outcome for this year
  normalized_data_year <- normalized_data %>%
    mutate(post_year = ifelse(entry_year_30 >= year, 1, 0))
  
  # Run RDD
  rd_result <- rdrobust(
    y = normalized_data_year$post_year,
    x = normalized_data_year$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd"
  )
  
  placebo_results <- rbind(placebo_results, data.frame(
    year = year,
    estimate = rd_result$coef[1],
    se = rd_result$se[1],
    pvalue = rd_result$pv[1],
    ci_lower = rd_result$ci[1,1],
    ci_upper = rd_result$ci[1,2]
  ))
  
  cat(sprintf("Year %d: Est = %.3f (SE = %.3f), p = %.4f %s\n",
              year, rd_result$coef[1], rd_result$se[1], rd_result$pv[1],
              ifelse(year == 2002, "← Main year", "")))
}

# Plot placebo results
p_placebo <- ggplot(placebo_results, aes(x = year, y = estimate)) +
  geom_line(color = "gray50") +
  geom_point(aes(color = factor(year == 2002)), size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, alpha = 0.5) +
  geom_vline(xintercept = 2002, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"), guide = FALSE) +
  labs(title = "Placebo Tests: RD Estimates for Different Year Thresholds",
       x = "Year Threshold",
       y = "RD Estimate",
       caption = "2002 (red) is the actual policy change year") +
  theme_minimal()

ggsave("output/robustness_placebo_national.pdf", p_placebo, width = 10, height = 6)

# ============================================================================
# TEST 4: ALTERNATIVE KERNELS
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 4: ALTERNATIVE KERNELS\n")
cat("------------------------------------\n\n")

kernels <- c("triangular", "uniform", "epanechnikov")
kernel_results <- data.frame()

for (kern in kernels) {
  rd_result <- rdrobust(
    y = normalized_data$post_2002,
    x = normalized_data$normalized_iml,
    c = 0,
    kernel = kern,
    bwselect = "mserd"
  )
  
  kernel_results <- rbind(kernel_results, data.frame(
    kernel = kern,
    estimate = rd_result$coef[1],
    se = rd_result$se[1],
    pvalue = rd_result$pv[1]
  ))
  
  cat(sprintf("%-12s: Est = %.3f (SE = %.3f), p = %.4f\n",
              kern, rd_result$coef[1], rd_result$se[1], rd_result$pv[1]))
}

# ============================================================================
# TEST 5: POLYNOMIAL ORDERS
# ============================================================================

cat("\n------------------------------------\n")
cat("TEST 5: POLYNOMIAL ORDERS\n")
cat("------------------------------------\n\n")

poly_orders <- 1:3
poly_results <- data.frame()

for (p_order in poly_orders) {
  rd_result <- rdrobust(
    y = normalized_data$post_2002,
    x = normalized_data$normalized_iml,
    c = 0,
    kernel = "triangular",
    bwselect = "mserd",
    p = p_order
  )
  
  poly_results <- rbind(poly_results, data.frame(
    order = p_order,
    estimate = rd_result$coef[1],
    se = rd_result$se[1],
    pvalue = rd_result$pv[1]
  ))
  
  cat(sprintf("Order %d: Est = %.3f (SE = %.3f), p = %.4f\n",
              p_order, rd_result$coef[1], rd_result$se[1], rd_result$pv[1]))
}

# ============================================================================
# SUMMARY TABLE
# ============================================================================

cat("\n========================================\n")
cat("ROBUSTNESS CHECK SUMMARY\n")
cat("========================================\n\n")

# Compile summary
summary_robust <- data.frame(
  Test = c("Main Result",
           "Donut 0.10",
           "Donut 0.20",
           "BW: CER",
           "Kernel: Uniform",
           "Polynomial: 2nd"),
  Estimate = c(main_result$coef[1],
               donut_results$estimate[donut_results$donut_size == 0.10],
               donut_results$estimate[donut_results$donut_size == 0.20],
               bw_results$estimate[bw_results$method == "cerrd"],
               kernel_results$estimate[kernel_results$kernel == "uniform"],
               poly_results$estimate[poly_results$order == 2]),
  SE = c(main_result$se[1],
         donut_results$se[donut_results$donut_size == 0.10],
         donut_results$se[donut_results$donut_size == 0.20],
         bw_results$se[bw_results$method == "cerrd"],
         kernel_results$se[kernel_results$kernel == "uniform"],
         poly_results$se[poly_results$order == 2])
)

summary_robust$Estimate <- round(summary_robust$Estimate, 3)
summary_robust$SE <- round(summary_robust$SE, 3)

print(summary_robust)

# Save all results
save(donut_results, bw_results, placebo_results, kernel_results, poly_results,
     file = "output/national_robustness_results.RData")

write.csv(summary_robust, "output/national_robustness_summary.csv", row.names = FALSE)

cat("\n✓ All robustness checks complete\n")
cat("✓ Results saved to output/national_robustness_results.RData\n")
cat("✓ Plots saved to output/robustness_*.pdf\n")
