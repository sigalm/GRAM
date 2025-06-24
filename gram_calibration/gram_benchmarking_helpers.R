### Benchmarking helper functions

benchmark_values <<- tribble(
  ~Variable,                    ~`Benchmark Value`,     ~Source,
  "n",                          "",                     "",
  "AGE",                        "",                     "",
  "SEX: Male",                  "49.0%",                "ACS data",
  "SEX: Female",                "51.0%",                "ACS data",
  "RACEETH: Non-Hispanic White","64.0%",                "ACS data",
  "RACEETH: Non-Hispanic Black","14.0%",                "ACS data",
  "RACEETH: Hispanic",          "22.0%",                "ACS data",
  "EDU",                        "13.6 (??)",            "UNDP-HDI",
  "INCOME: Low",                "5.0%",                 "ACS data",
  "INCOME: Medium",             "17.0%",                "ACS data",
  "INCOME: High",               "78.0%",                "ACS data",
  "MEDBUR",                     "2.25 (??)",            "Mossadeghi et al 2023",
  "APOE4: Not Carrier",         "75.0%",                "Di Battista et al 2016",
  "APOE4: Carrier",             "25.0%",                "Di Battista et al 2016",
  "SYN: Cognitively Healthy",   "",                     ""
)


lifetable <- read.csv("gram_data/mortality/table_prob_die_next_year.csv")
lifetable <<- lifetable %>%
  mutate(rate = - log(1 - qx))

benchmark_prev_by_age <<- tribble(
  ~ age,      ~ raceeth,      ~ condition,    ~ source,     ~ prev,    ~ ci_lo,    ~ ci_hi,
  67,         "Overall",      "dem",          "Manly",      0.03,       0.01,       0.04,
  72,         "Overall",      "dem",          "Manly",      0.04,       0.02,       0.06,
  77,         "Overall",      "dem",          "Manly",      0.09,       0.06,       0.11,
  82,         "Overall",      "dem",          "Manly",      0.18,       0.14,       0.22,
  87,         "Overall",      "dem",          "Manly",      0.26,       0.20,       0.31,
  95,         "Overall",      "dem",          "Manly",      0.35,       0.28,       0.43,
  
  # 67,         "Overall",    "mci",          "Manly",       0.22,       0.18,       0.25,
  # 72,         "Overall",    "mci",          "Manly",       0.20,       0.17,       0.24,
  # 77,         "Overall",    "mci",          "Manly",       0.21,       0.18,       0.24,
  # 82,         "Overall",    "mci",          "Manly",       0.25,       0.21,       0.29,
  # 87,         "Overall",    "mci",          "Manly",       0.22,       0.17,       0.27,
  # 95,         "Overall",    "mci",          "Manly",       0.27,       0.20,       0.35,
  
  65,         "Overall",      "mci",          "Bai",         0.115,      0.076,      0.161,
  75,         "Overall",      "mci",          "Bai",         0.158,      0.118,      0.202,
  88,         "Overall",      "mci",          "Bai",         0.213,      0.163,      0.268
)


proportion_instit <<- tribble(
  ~ age,    ~proportion,   ~prev_instit,
  67,       0.013,         0.40,
  72,       0.013,         0.45,
  77,       0.10,          0.50,
  82,       0.10,          0.55,
  87,       0.30,          0.60,
  94,       0.30,          0.70
)

benchmark_prev_by_age_instit <<- benchmark_prev_by_age %>%
  left_join(proportion_instit, by = "age") %>%
  mutate(prev_adj = proportion * prev_instit + (1-proportion) * prev)

flextable(benchmark_prev_by_age_instit)

benchmark_prev_by_race <<- data.frame(
  age = rep(74, times = 3),
  raceeth = c("NHW","NHB","Hisp"),
  condition = c(rep("dem_manly", times = 3), rep("mci_manly", times = 3)),
  prev = c(0.11, 0.15, 0.10,
           0.23, 0.22, 0.28)
)

benchmark_reside_time <<- data.frame(
  age_group = "Benchmark",
  condition = c("benchmark_mci", "benchmark_dem"),
  duration = c(3.85, 7)
)

# Note that the benchmark dementia duration is a guess based on conversations with Kelly and Kate.

benchmark_age_of_onset <<- 71.5



compare_mortality <- function(sim, description, n) {
  mortality <- as.data.frame(sim$aggregated_results_totpop$state_trace) %>%
    mutate(new_dth = dth - lag(dth),  # new deaths per cycle
           denom = 1 - lag(dth),       # those alive at the start of cycle
           model_rate = new_dth / denom) %>%
    mutate(model_rate = replace_na(model_rate, 0))
  
  mort_compare <- cbind(lifetable, mortality) %>%
    mutate(benchmark_rate_1000 = rate * 1000,
           model_rate_1000 = model_rate * 1000,
           residual = model_rate_1000 - benchmark_rate_1000) %>%
    mutate(
      cum_model_rate = mortality[ ,"dth"] * 1000,
      cum_benchmark_rate = (1 - cumprod(1-qx)) * 1000,
      cum_residual = cum_model_rate - cum_benchmark_rate
    )
  
  fig_mort_compare <- ggplot(mort_compare, aes(x = age)) +
    geom_line(aes(y = benchmark_rate_1000, color = "Benchmark", linetype = "Benchmark"), linewidth = 1) +
    geom_line(aes(y = model_rate_1000, color = "Model", linetype = "Model"), linewidth = 1) +
    scale_color_manual(values = c("Benchmark" = "grey", "Model" = "black")) +
    scale_linetype_manual(values = c("Benchmark" = "solid", "Model" = "dashed")) +
    labs(title = "Mortality Rate",
         x = "Age", y = "Mortality Rate", color = "Source", linetype = "Source",
         caption = "Mortality Rate per 1000 Individuals") +
    theme_minimal(base_size = 16)
  
  fig_mort_compare_residuals <- ggplot(mort_compare, aes(x = age, y = residual)) +
    geom_line() + 
    geom_point() + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
    labs(title = "Residuals",
         x = "Age",
         y = "Residual") +
    theme_minimal(base_size = 16)
  
  fig_cum_mort_compare <- ggplot(mort_compare, aes(x = age)) +
    geom_line(aes(y = cum_benchmark_rate, color = "Benchmark", linetype = "Benchmark"), size = 1) +
    geom_line(aes(y = cum_model_rate, color = "Model", linetype = "Model"), size = 1) +
    scale_color_manual(values = c("Benchmark" = "grey", "Model" = "black")) +
    scale_linetype_manual(values = c("Benchmark" = "solid", "Model" = "dashed")) +
    labs(title = "Cumulative Mortality Rate",
         x = "Age", y = "Mortality Rate", color = "Source", linetype = "Source",
         caption = "Mortality Rate per 1000 Individuals") +
    theme_minimal(base_size = 16)
  
  fig_cum_mort_compare_residuals <- ggplot(mort_compare, aes(x = age, y = cum_residual)) +  
    geom_line() + 
    geom_point() + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
    labs(title = "Residuals",
         x = "Age",
         y = "Residual") +
    theme_minimal(base_size = 16)
  
  result_plot <- (fig_mort_compare + fig_mort_compare_residuals) / (fig_cum_mort_compare + fig_cum_mort_compare_residuals) +
    plot_annotation(title = "Modeled vs. Empirical Mortality Rate", subtitle = paste0(description, "\nN = ", n))
  
  return(invisible(list(result_plot = result_plot,
                        dat = mort_compare,
                        annual_mort = list(fig_mort_compare, fig_mort_compare_residuals),
                        cumulative_mort = list(fig_cum_mort_compare, fig_cum_mort_compare_residuals)
  )))
}



 
# Stratify prevalence by age and any variable
stratify_prevalence_by <- function(sim, strat_var, strat_labels = NULL, strat_cutoffs = NULL) {
  
  # Extract stratification variable
  strat_vec <- as.vector(sim$output[, strat_var,])
  # Bin if cutoffs provided
  if (!is.null(strat_cutoffs)) {
    strat_binned <- cut(strat_vec, breaks = strat_cutoffs, include.lowest = TRUE, right = FALSE)
    strat <- strat_binned
  } else {
    strat <- strat_vec
  }

  # Extract variables
  df <- data.frame(
    id = rep(1:dim(sim$output)[3], times = dim(sim$output)[1]),
    age = as.vector(sim$output[,"AGE",]),
    alive = as.vector(sim$output[,"ALIVE",]),
    sev = as.vector(sim$output[,"SEV",]),
    strat = strat
  ) %>%
    filter(alive == 1) %>%
    mutate(condition = case_when(
      sev == 0 ~ "MCI",
      sev %in% c(1,2,3) ~ "Dementia",
      TRUE ~ "Healthy"
    ))

  # Apply labels if provided (for factor or binned stratification)
  if (!is.null(strat_labels)) {
    if (!is.null(strat_cutoffs)) {
      df$strat <- factor(df$strat, labels = strat_labels)
    } else {
      df$strat <- factor(df$strat, levels = seq_along(strat_labels) - 1, labels = strat_labels)
    }
  }

  prev_strat <- df %>%
    group_by(age, strat) %>%
    summarise(
      tot_alive = n(),
      mci = sum(condition == "MCI") / tot_alive,
      dem = sum(condition == "Dementia") / tot_alive,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mci, dem), names_to = "condition", values_to = "prev")

  plot_prev <- ggplot(prev_strat, aes(x = age, y = prev, color = condition)) +
    geom_line() +
    facet_wrap(~ strat) +
    labs(title = paste0("Prevalence by ", strat_var), x = "Age", y = "Prevalence", color = "Condition") +
    scale_color_manual(
      labels = c("dem" = "Dementia", "mci" = "MCI"),
      values = c("dem" = "darkgreen", "mci" = "hotpink")) +
    theme_minimal(base_size = 14)

  return(invisible(list(plot_prev = plot_prev,
                        dat = prev_strat)))
}

compare_prevalence <- function(sim, description, n) {
  df <- data.frame(
    id = rep(1:dim(sim$output)[3], times = dim(sim$output)[1]),
    age = as.vector(sim$output[,"AGE",]),
    alive = as.vector(sim$output[,"ALIVE",]),
    sev = as.vector(sim$output[,"SEV",])
  ) %>%
    filter(alive == 1) %>%
    mutate(condition = case_when(
      sev == 0 ~ "MCI",
      sev %in% c(1,2,3) ~ "Dementia",
      TRUE ~ "Healthy"
    ))
  
  prev_overall <- df %>%
    group_by(age) %>%
    summarise(
      tot_alive = n(),
      mci = sum(condition == "MCI") / tot_alive,
      dem = sum(condition == "Dementia") / tot_alive,
      .groups = "drop"
    ) %>%
    pivot_longer(cols = c(mci, dem), names_to = "condition", values_to = "prev") %>%
    mutate(cond = case_when(
      condition %in% c("dem", "Dementia") ~ "Dementia",
      condition %in% c("mci", "MCI") ~ "MCI",
      TRUE ~ as.character(condition)),
      source = "Model"
    )
  
  # For plotting, harmonize 'condition' to "Dementia"/"MCI" for color
  
  bench_prev <- benchmark_prev_by_age %>% mutate(
    cond = case_when(
      condition == "dem" ~ "Dementia",
      condition == "mci" ~ "MCI",
      TRUE ~ as.character(condition)
    ),
    source = as.character(source)
  )
  
  # Combine for plotting
  plot_prev <- bind_rows(
    prev_overall %>% select(age, prev, source, cond),
    bench_prev %>% select(age, prev, source, cond, ci_lo, ci_hi)
  )
  
  # Plot
  # Unify benchmark sources for legend
  plot_prev <- plot_prev %>% mutate(source_for_legend = ifelse(source == "Model", "Model", "Benchmark"))
  
  fig_prev_by_age <- ggplot() +
    geom_line(data = plot_prev %>% filter(source_for_legend == "Model"),
              aes(x = age, y = prev, color = cond), size = 1) +
    geom_point(data = plot_prev %>% filter(source_for_legend == "Benchmark"),
               aes(x = age, y = prev, color = cond), size = 2) +
    geom_errorbar(data = plot_prev %>% filter(source_for_legend == "Benchmark" & !is.na(ci_lo)),
                  aes(x = age, ymin = ci_lo, ymax = ci_hi, color = cond), width = 0.5) +
    scale_color_manual(
      name = NULL,
      values = c("Dementia" = "darkgreen", "MCI" = "hotpink")) +
    labs(title = "Prevalence of Cognitive Impairment by Age",
         subtitle = paste0(description, "\nN = ", n),
         x = "Age",
         y = "Prevalence") +
    ylim(0, 1) +
    theme_minimal(base_size = 16)
  
  
  return(invisible(list(fig_prev_by_age = fig_prev_by_age,
                        dat = prev_overall)))
}


compare_reside_time <- function(sim, description, n) {
  reside_time <- as.data.frame(sim$aggregated_results_totpop$reside_time$noncensored) %>%
    # mutate(dem = mil + mod + sev) %>%
    select(-mil, -mod, - sev) %>%
    pivot_longer(cols = c(mci, any_dem), names_to = "condition", values_to = "duration")
  
  age_bins <- cut(sim$aggregated_results_totpop$age_at_onset, 
                  breaks = seq(50, 100, by = 5), 
                  right = FALSE, 
                  include.lowest = TRUE)
  
  age_group_weights <- as.data.frame(table(age_bins)) %>%
    rename(age_group = age_bins, count = Freq) %>%
    mutate(weight = count / sum(count),
           age_group = unique(reside_time$age_group)[-11])
  
  reside_time_adjusted <- reside_time %>%
    filter(age_group != "Overall") %>%
    left_join(age_group_weights, by = "age_group") %>%
    mutate(weighted_duration = duration * weight) %>%
    group_by(condition) %>%
    summarise(duration = sum(weighted_duration, na.rm = TRUE)) %>%
    mutate(age_group = "Overall_adj") %>%
    select(age_group, condition, duration)
  
  reside_time <- reside_time %>%
    rbind(reside_time_adjusted)
  
  # Remove rows where age_group == 'Overall'
  reside_time <- reside_time[reside_time$age_group != "Overall", ]
  # Rename 'Overall_adj' to 'Overall'
  reside_time$age_group[reside_time$age_group == "Overall_adj"] <- "Overall"
  
  fig_reside_time <- ggplot() +
    geom_col(data = reside_time[reside_time$age_group != "Overall", ],
             aes(x = age_group, y = duration, fill = condition), width = 0.6) +
    geom_col(data = reside_time[reside_time$age_group == "Overall", ],
             aes(x = age_group, y = duration, fill = condition), width = 0.6, alpha = 0.85) +
    geom_col(data = benchmark_reside_time, aes(x = age_group, y = duration, fill = condition), width = 0.6) +
    labs(title = "Reside Time in MCI and Dementia by Age of Onset",
         subtitle = paste0("N = ", n),
         x = "Age of Onset",
         y = "Duration (years)",
         fill = NULL) +
    scale_fill_manual(
      labels = c("any_dem" = "Model: Dementia", "mci" = "Model: MCI", 
                 "benchmark_dem" = "Benchmark: Dementia", "benchmark_mci" = "Benchmark: MCI"),
      values = c("any_dem" = "darkgreen", "mci" = "hotpink", "benchmark_dem" = "lightgreen", "benchmark_mci" = "pink")) +
    theme_minimal(base_size = 16) +
    theme(axis.text.x = element_text(angle = 30, hjust = 0.5))
  
  return(invisible(list(fig_reside_time = fig_reside_time,
                        dat = reside_time)))
  
}


compare_age_onset <- function(sim, description, n) {
  
  age_onset <- sim$aggregated_results_totpop$age_at_onset
  avg_age_onset <- mean(age_onset, na.rm = TRUE)
  df_age_onset <- data.frame(age_onset = age_onset, apoe4 = sim$output[1,"APOE4",],
                             raceeth = sim$output[1,"RACEETH",])
  
  
  age_onset_compare <- data.frame(
    type = c("Model", "Benchmark"),
    value = c(avg_age_onset, benchmark_age_of_onset)
  )
  
  fig_age_onset <- ggplot(data = age_onset_compare, aes(x = type, y = value)) +
    geom_col(aes(fill = type), width = 0.6) + 
    geom_text(aes(label = round(value, 1)), vjust = -0.5) + 
    labs(title = "Modeled vs. Benchmark Mean Age of MCI Onset",
         subtitle = paste0(description, "\nN = ", n),
         x = NULL,
         y = "Age") + 
    scale_fill_manual(values = c("Benchmark" = "grey", "Model" = "darkgreen")) +
    expand_limits(y = max(age_onset_compare$value) * 1.1) +
    guides(fill = "none") +
    theme_minimal()
  
  avg_onset_by_apoe4 <- df_age_onset %>%
    filter(!is.na(age_onset)) %>%
    group_by(apoe4) %>%
    summarize(mean = mean(age_onset))
  avg_onset_by_raceeth <- df_age_onset %>%
    filter(!is.na(age_onset)) %>%
    group_by(raceeth) %>%
    summarize(mean = mean(age_onset))
  
  
  return(invisible(list(fig_age_onset = fig_age_onset,
                        dat = age_onset_compare)))
}

run_benchmarking <- function(l.inputs, description = NULL, sample = NULL) {
  sim <- f.wrap_run(l.inputs = l.inputs, microdata = sample)
  sim_desc <- description
  n <- l.inputs[["n.ind"]]
  
  mort <- compare_mortality(sim, sim_desc, n)
  prev <- compare_prevalence(sim, sim_desc, n)
  reside_time <- compare_reside_time(sim, sim_desc, n)
  age_onset <- compare_age_onset(sim, sim_desc, n)
  
  return(invisible(list(sim = sim, mort = mort, prev = prev, reside_time = reside_time, age_onset = age_onset)))
}