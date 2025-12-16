######################################## GRAM HELPER FUNCTIONS: OUTPUT FORMATTERS ########################################

f.format_reside_time_table <- function(reside_time_data, scenario_title = NULL) {
  table <- flextable(reside_time_data) %>%
    set_header_labels(
      age_group = "Age at MCI Onset",
      mci = "MCI",
      mil = "Mild Dementia",
      mod = "Moderate Dementia",
      sev = "Severe Dementia"
    ) %>%
    add_header_row(
      values = c("", "Average Reside Time in Years, by Age at Onset"),
      colwidths = c(1, 4)
    ) %>%
    add_header_lines(values = scenario_title) %>%
    bold(j = 1) %>%
    theme_vanilla() %>%
    autofit() %>%
    align(j = 2:5, align = "center", part = "all") %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all")
  
  return(table) 
}

f.format_prevalence_table <- function(prevalence_data, denom = 1000, scenario_title = NULL) {
  
  dat <- prevalence_data
  dat[ ,-1] <- round(dat[ ,-1] * denom, digits = 2)
  
  table <- flextable(dat) %>%
    set_header_labels(
      age_group = "Age",
      avg_mci = "MCI",
      avg_mil = "Mild Dementia",
      avg_mod = "Moderate Dementia",
      avg_sev = "Severe Dementia"
    ) %>%
    add_header_row(
      values = c("", paste0("Prevalence by Severity by Age, per ", denom)),
      colwidths = c(1, 4)
    ) %>%
    add_header_lines(values = scenario_title) %>%
    bold(j = 1) %>%
    theme_vanilla() %>%
    autofit() %>%
    align(j = 2:5, align = "center", part = "all") %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all")
  
}