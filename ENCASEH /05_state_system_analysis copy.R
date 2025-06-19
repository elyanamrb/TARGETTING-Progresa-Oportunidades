# 04_state_system_analysis.R
# Analysis of scoring systems by state (ENT_FED)

library(tidyverse)
library(ggplot2)
library(scales)

# Load prepared data
prepared_data <- readRDS("prepared_encaseh_data.rds")

# Create directory for state analysis
dir.create("state_analysis", showWarnings = FALSE)

# Mexican state codes and names (for reference)
state_names <- data.frame(
  ENT_FED = 1:32,
  State_Name = c(
    "Aguascalientes", "Baja California", "Baja California Sur", "Campeche",
    "Coahuila", "Colima", "Chiapas", "Chihuahua", "Ciudad de México",
    "Durango", "Guanajuato", "Guerrero", "Hidalgo", "Jalisco",
    "México", "Michoacán", "Morelos", "Nayarit", "Nuevo León",
    "Oaxaca", "Puebla", "Querétaro", "Quintana Roo", "San Luis Potosí",
    "Sinaloa", "Sonora", "Tabasco", "Tamaulipas", "Tlaxcala",
    "Veracruz", "Yucatán", "Zacatecas"
  )
)

# Function to analyze scoring systems by state
analyze_systems_by_state <- function(data, year_label) {
  cat(sprintf("\n=== State-Level Scoring System Analysis for Year %s ===\n", year_label))
  
  # Add scoring system classification
  data <- data %>%
    mutate(
      scoring_system = ifelse(PUN_INC > 100, "System 2 (High)", "System 1 (Normal)"),
      ENT_FED = as.numeric(as.character(ENT_FED))
    )
  
  # Calculate statistics by state
  state_summary <- data %>%
    group_by(ENT_FED) %>%
    summarise(
      total_households = n(),
      n_system1 = sum(scoring_system == "System 1 (Normal)"),
      n_system2 = sum(scoring_system == "System 2 (High)"),
      pct_system2 = n_system2 / total_households * 100,
      mean_pun_inc = mean(PUN_INC, na.rm = TRUE),
      median_pun_inc = median(PUN_INC, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(state_names, by = "ENT_FED") %>%
    arrange(desc(pct_system2))
  
  # Identify states predominantly using System 2
  system2_states <- state_summary %>%
    filter(pct_system2 > 50) %>%
    pull(State_Name)
  
  if (length(system2_states) > 0) {
    cat(sprintf("\nStates predominantly using System 2 (>50%%):\n"))
    for (state in system2_states) {
      pct <- state_summary$pct_system2[state_summary$State_Name == state]
      cat(sprintf("- %s: %.1f%%\n", state, pct))
    }
  } else {
    cat("\nNo states predominantly using System 2\n")
  }
  
  # States with mixed systems (10-90%)
  mixed_states <- state_summary %>%
    filter(pct_system2 >= 10 & pct_system2 <= 90) %>%
    arrange(desc(pct_system2))
  
  if (nrow(mixed_states) > 0) {
    cat(sprintf("\nStates with mixed systems (10-90%% System 2):\n"))
    for (i in 1:min(10, nrow(mixed_states))) {
      cat(sprintf("- %s: %.1f%%\n", mixed_states$State_Name[i], mixed_states$pct_system2[i]))
    }
  }
  
  # Create visualizations
  # 1. Bar chart of System 2 usage by state
  state_plot <- ggplot(state_summary, aes(x = reorder(State_Name, pct_system2), y = pct_system2)) +
    geom_bar(stat = "identity", aes(fill = pct_system2)) +
    scale_fill_gradient2(low = "lightblue", mid = "orange", high = "red", 
                         midpoint = 50, limits = c(0, 100)) +
    coord_flip() +
    geom_hline(yintercept = 50, linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = sprintf("Percentage of Households Using System 2 by State - %s", year_label),
         subtitle = "System 2 = High scores (>100), System 1 = Normal scores (≤100)",
         x = "State",
         y = "Percentage Using System 2",
         fill = "% System 2") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # 2. Map-style visualization (if we have enough variation)
  if (sum(state_summary$pct_system2 > 0) > 1) {
    # Create a heatmap-style plot
    map_data <- state_summary %>%
      mutate(
        x = (ENT_FED - 1) %% 8,
        y = 7 - floor((ENT_FED - 1) / 8)
      )
    
    heatmap_plot <- ggplot(map_data, aes(x = x, y = y)) +
      geom_tile(aes(fill = pct_system2), color = "white", size = 0.5) +
      geom_text(aes(label = sprintf("%s\n%.0f%%", 
                                    substr(State_Name, 1, 8), 
                                    pct_system2)), 
                size = 2.5) +
      scale_fill_gradient2(low = "lightblue", mid = "orange", high = "red", 
                           midpoint = 50, limits = c(0, 100)) +
      labs(title = sprintf("State Grid: System 2 Usage - %s", year_label),
           fill = "% System 2") +
      theme_void() +
      theme(legend.position = "bottom")
  } else {
    heatmap_plot <- NULL
  }
  
  # 3. Scatter plot of mean PUN_INC vs % System 2
  scatter_plot <- ggplot(state_summary, aes(x = pct_system2, y = mean_pun_inc)) +
    geom_point(aes(size = total_households), alpha = 0.6, color = "darkblue") +
    geom_text(data = filter(state_summary, pct_system2 > 20 | mean_pun_inc > 50),
              aes(label = substr(State_Name, 1, 10)), 
              vjust = -0.5, hjust = 0.5, size = 3) +
    scale_size_continuous(range = c(2, 10), labels = scales::comma) +
    labs(title = sprintf("State Characteristics by Scoring System - %s", year_label),
         x = "Percentage Using System 2",
         y = "Mean PUN_INC",
         size = "Total\nHouseholds") +
    theme_minimal()
  
  # Save plots
  ggsave(sprintf("state_analysis/state_system2_usage_%s.pdf", year_label), 
         state_plot, width = 10, height = 12)
  
  if (!is.null(heatmap_plot)) {
    ggsave(sprintf("state_analysis/state_heatmap_%s.pdf", year_label), 
           heatmap_plot, width = 10, height = 8)
  }
  
  ggsave(sprintf("state_analysis/state_scatter_%s.pdf", year_label), 
         scatter_plot, width = 10, height = 8)
  
  # Save detailed table
  write.csv(state_summary, 
            sprintf("state_analysis/state_summary_%s.csv", year_label), 
            row.names = FALSE)
  
  return(list(
    state_summary = state_summary,
    system2_states = system2_states,
    plots = list(state_plot = state_plot, 
                 heatmap_plot = heatmap_plot,
                 scatter_plot = scatter_plot)
  ))
}

# Function to compare state patterns across years
compare_states_across_years <- function(all_state_results) {
  # Extract data for comparison
  comparison_data <- data.frame()
  
  for (year_name in names(all_state_results)) {
    year <- gsub("y", "", year_name)
    state_data <- all_state_results[[year_name]]$state_summary %>%
      select(ENT_FED, State_Name, pct_system2) %>%
      mutate(Year = as.numeric(year))
    
    comparison_data <- rbind(comparison_data, state_data)
  }
  
  # Create evolution plot for states with significant System 2 usage
  significant_states <- comparison_data %>%
    group_by(State_Name) %>%
    summarise(max_system2 = max(pct_system2)) %>%
    filter(max_system2 > 10) %>%
    pull(State_Name)
  
  if (length(significant_states) > 0) {
    evolution_plot <- comparison_data %>%
      filter(State_Name %in% significant_states) %>%
      ggplot(aes(x = Year, y = pct_system2, color = State_Name, group = State_Name)) +
      geom_line(size = 1) +
      geom_point(size = 2) +
      labs(title = "Evolution of System 2 Usage by State (2000-2004)",
           subtitle = "Only showing states that had >10% System 2 usage at some point",
           x = "Year",
           y = "Percentage Using System 2") +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggsave("state_analysis/state_system2_evolution.pdf", 
           evolution_plot, width = 12, height = 8)
  }
  
  # Create summary table of transitions
  transition_summary <- comparison_data %>%
    filter(Year %in% c(2000, 2001)) %>%
    select(State_Name, Year, pct_system2) %>%
    pivot_wider(names_from = Year, values_from = pct_system2, names_prefix = "Year_") %>%
    mutate(
      Change_2000_2001 = Year_2001 - Year_2000,
      Status_2000 = case_when(
        Year_2000 > 90 ~ "Predominantly System 2",
        Year_2000 > 10 ~ "Mixed Systems",
        Year_2000 > 0 ~ "Minimal System 2",
        TRUE ~ "System 1 Only"
      ),
      Status_2001 = case_when(
        Year_2001 > 90 ~ "Predominantly System 2",
        Year_2001 > 10 ~ "Mixed Systems",
        Year_2001 > 0 ~ "Minimal System 2",
        TRUE ~ "System 1 Only"
      )
    ) %>%
    arrange(desc(Year_2001))
  
  write.csv(transition_summary, 
            "state_analysis/state_system_transitions_2000_2001.csv", 
            row.names = FALSE)
  
  return(transition_summary)
}

# Run analysis for each year
all_state_results <- list()

for (year_name in c("y2000", "y2001", "y2002", "y2004")) {
  year_label <- gsub("y", "", year_name)
  all_state_results[[year_name]] <- analyze_systems_by_state(
    prepared_data[[year_name]], 
    year_label
  )
}

# Compare across years
transition_summary <- compare_states_across_years(all_state_results)

# Print key findings
cat("\n\n=== KEY FINDINGS: State-Level Scoring Systems ===\n")

# Year 2001 focus
cat("\n=== Year 2001 - Critical Transition Year ===\n")
year_2001_summary <- all_state_results$y2001$state_summary %>%
  filter(pct_system2 > 0) %>%
  arrange(desc(pct_system2))

cat("\nStates still using System 2 in 2001:\n")
for (i in 1:min(15, nrow(year_2001_summary))) {
  cat(sprintf("%2d. %-20s: %6.1f%% (n=%d)\n", 
              i,
              year_2001_summary$State_Name[i], 
              year_2001_summary$pct_system2[i],
              year_2001_summary$n_system2[i]))
}

# Transition patterns
cat("\n\nMajor transitions from 2000 to 2001:\n")
major_transitions <- transition_summary %>%
  filter(abs(Change_2000_2001) > 20) %>%
  arrange(Change_2000_2001)

if (nrow(major_transitions) > 0) {
  for (i in 1:nrow(major_transitions)) {
    cat(sprintf("- %s: %.1f%% → %.1f%% (change: %+.1f%%)\n",
                major_transitions$State_Name[i],
                major_transitions$Year_2000[i],
                major_transitions$Year_2001[i],
                major_transitions$Change_2000_2001[i]))
  }
}

cat("\nAnalysis complete! Check 'state_analysis/' directory for detailed results.\n")