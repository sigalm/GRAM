# ********************************************************************************
# ======= GRAM SCENARIOS FOR JIM SOUTH AFRICA PRESENTATION - 2024 NOVEMBER ======= 
# ********************************************************************************
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")


l.inputs[["r.CDRslow_mean"]] <- 0.3
l.inputs1 <- l.inputs


# Scenario 1 - Natural progression, no treatment
l.inputs1[["scenario"]] <- "Natural history, no intervention"
scenario1 <- f.wrap_run(l.inputs1)
fig_scenario1 <- scenario1$fig.progression$fig.progression_true
reside_time_scenario1 <- format_reside_timefig.progression_truereside_time_scenario1 <- format_reside_time_table(scenario1$aggregated_results_totpop$reside_time_by_age_of_onset)
prevalence_scenario1 <- format_prevalence_table(scenario1$aggregated_results_totpop$prevalence_grouped)

# Scenario 2 - Perfect information: DMT introduced for MCI for max 3 years or moderate dementia
l.inputs2 <- l.inputs1
l.inputs2[["scenario"]] <- "Start DMT at MCI"
l.inputs2[["Tx"]] <- 1
scenario2 <-  f.wrap_run(l.inputs2)
fig_scenario2 <- scenario2$fig.progression
reside_time_scenario2 <- format_reside_time_table(scenario2$aggregated_results_totpop$reside_time_by_age_of_onset)
prevalence_scenario2 <- format_prevalence_table(scenario2$aggregated_results_totpop$prevalence_grouped)


# Scenario 3 - Imperfect information: DMT introduced for mild dementia for max 3 years of moderate dementia
l.inputs3 <- l.inputs1
l.inputs3[["scenario"]] <- "Start DMT at mild dementia"
l.inputs3[["Tx"]] <- 1
l.inputs3[["sens_BHA"]][1] <- 0
l.inputs3[["spec_BHA"]] <- 1
scenario3 <- f.wrap_run(l.inputs3)
fig_scenario3 <- scenario3$fig.progression
reside_time_scenario3 <- format_reside_time_table(scenario3$aggregated_results_totpop$reside_time_by_age_of_onset)
prevalence_scenario3 <- format_prevalence_table(scenario3$aggregated_results_totpop$prevalence_grouped)


# Scenario 4 - Hypothetical prevention intervention that reduces the probability of developing MCI
l.inputs4 <- l.inputs1
l.inputs4[["scenario"]] <- "Hypothetical MCI prevention intervention with 30% effectiveness, no DMT"
l.inputs4[["rr.Px_mci"]] <- 1-0.3
scenario4 <- f.wrap_run(l.inputs4)
fig_scenario4 <- scenario4$fig.progression
reside_time_scenario4 <- format_reside_time_table(scenario4$aggregated_results_totpop$reside_time_by_age_of_onset)
prevalence_scenario4 <- format_prevalence_table(scenario4$aggregated_results_totpop$prevalence_grouped)


# Scenario 5 - Everyone starts with MCI and at age 70
l.inputs5 <- l.inputs
l.inputs5[["scenario"]] <- "everyone starts with MCI"
l.inputs5[["AGE_start_mean"]] <- 70
l.inputs5[["n.cycle"]] <- 30
l.inputs5[["p.SYN_start"]] <- c(0,1)
l.inputs5[["r.CDRslow_mean"]] <- 0.3 # original from lit, we then reduced to 0.3
scenario5 <- f.wrap_run(l.inputs5)
fig_scenario5 <- scenario5$fig.progression$fig.progression_true
reside_time_scenario5 <- format_reside_time_table(scenario5$aggregated_results_totpop$reside_time_by_age_of_onset)



# Save plots as png files
ggsave("gram_results/2024_Nov_presentation/fig1_gram_natural_history.png", plot = fig_scenario1, width = 10, height = 6, dpi = 300)
ggsave("gram_results/2024_Nov_presentation/fig2_gram_dmt_at_mci.png", plot = fig_scenario2, width = 10, height = 6, dpi = 300)
ggsave("gram_results/2024_Nov_presentation/fig3_gram_dmt_at_milddem.png", plot = fig_scenario3, width = 10, height = 6, dpi = 300)
ggsave("gram_results/2024_Nov_presentation/fig4_gram_prevent_mci.png", plot = fig_scenario4, width = 10, height = 6, dpi = 300)

# Save tables to word
save_as_docx(
  "Natural history, no intervention" = reside_time_scenario1, 
  "Start DMT at MCI" = reside_time_scenario2,
  "Start DMT at mild dementia" = reside_time_scenario3,
  "MCI prevention" = reside_time_scenario4,
  path = "gram_results/2024_Nov_presentation/reside_time_tables.docx")
save_as_docx(
  "Natural history, no intervention" = prevalence_scenario1, 
  "Start DMT at MCI" = prevalence_scenario2,
  "Start DMT at mild dementia" = prevalence_scenario3,
  "MCI prevention" = prevalence_scenario4,
  path = "gram_results/2024_Nov_presentation/prevalence_tables.docx"
)
