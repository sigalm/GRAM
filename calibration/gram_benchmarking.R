---
title: "GRAM Simulation"
subtitle: "Calibration and Benchmarking"
author: "Sigal Maya"
date: "02/26/2025"
output: html_notebook
knit_root_dir: "/Users/smaya/Documents/GitHub/GRAM"
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = FALSE, include = FALSE, warning = FALSE, message = FALSE)
library(tidyverse)
getwd()
file.exists("gram_model/gram_01setup.r")
list.files("gram_model")
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
library(knitr)
```
```

# Run base model:
# 1. no considering of cognitive concerns for BHA testing
# 2. strict testing criteria

l.inputs1 <- l.inputs
l.inputs1[["n.cycle"]] <- 51
sim1 <- f.wrap_run(l.inputs1)

tab_prev_age <- f.format_prevalence_table(sim1$aggregated_results_totpop$prevalence_by_age)
tab_duration <- f.format_reside_time_table(sim1$aggregated_results_totpop$reside_time$noncensored)
avg_age_onset <- mean(sim1$aggregated_results_totpop$age_at_mci, na.rm = TRUE)

# Compare age-specific mortality to CDC-wonder all-cause mortality
mortality <- as.data.frame(sim1$aggregated_results_totpop$state_trace) %>%
  mutate(denom = rowSums(.) - lag(dth),
         new_dth = dth - lag(dth)) %>%
  mutate(mort = new_dth / denom) %>%
  select(mort)

# Re-generate `allcause_clean` dataframe from script "gram_data/calculate_non_dementia_death_rates.R"

mort_compare <- cbind(mortality, allcause_clean) %>%
  mutate(benchmark_mort = rate) %>%
  rename(model_mort = mort) %>%
  mutate(delta_abs = benchmark_mort - model_mort,
         delta_pct = (benchmark_mort - model_mort) / benchmark_mort) %>%
  select(age, model_mort, benchmark_mort, delta_abs, delta_pct)
max(abs(mort_compare$delta_pct), na.rm = TRUE)
min(abs(mort_compare$delta_pct), na.rm = TRUE)

kable(mort_compare)


# Review prevalence by race/eth
sim1$aggregated_results_totpop$prevalence_by_raceeth[25,,]  # at age 75, cycle 25
