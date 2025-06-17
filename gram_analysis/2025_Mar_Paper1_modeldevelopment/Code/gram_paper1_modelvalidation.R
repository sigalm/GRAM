# ********************************************************************************
# ======= GRAM PAPER 1: GRAM MODEL VALIDATION AND CALIBRATION ======= 
# ********************************************************************************

# TECHNICAL ##### 
# run only once, do not delete
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
source("gram_calibration/gram_benchmarking_helpers.r")
library(tableone)
sample1 <- readRDS("gram_data/acs_data/acs_sample_1.rds")
# sample2 <- readRDS("gram_data/acs_data/acs_sample_2.rds")
# sample3 <- readRDS("gram_data/acs_data/acs_sample_3.rds")
# ***************


# CALIBRATION ####
l.inputs1 <- l.inputs
l.inputs1[["n.cycle"]] <- 51
# l.inputs1[["m.lifetable"]] <- l.inputs[["m.lifetable"]] * 0.8
# 
# l.inputs1[["hr.mort_mci"]] <- 2
# l.inputs1[["hr.mort_mil"]] <- 3.5
# l.inputs1[["hr.mort_mod"]] <- 4.1
# l.inputs1[["hr.mort_sev"]] <- 10

# l.inputs1[["hr.mort_mil_age"]] <- l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(2, 2, 2) 
l.inputs1[["r.CDRslow_mean"]] <-  (seq(0, 1, length.out = 51)^1.5) * (1.5 *l.inputs[["r.CDRslow_mean"]])
# l.inputs1[["r.CDRfast_mean"]] <-  (seq(0.5, 1.5, length.out = 51)^2) * (2 * l.inputs[["r.CDRfast_mean"]])
l.inputs1[["p.HCARE_start"]] <- c(0.25,0.75)
l.inputs1[["n.ind"]] <- 50000
l.inputs1[["hr.mort_mci_age"]] <- c(1, 1, 1)
# l.inputs1[["hr.mort_mil_age"]] <- l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(2, 2, 2) 
l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(1, 1, 1) 
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 2.2 /4

mean(l.inputs1[["r.CDRslow_mean"]])
mean(l.inputs1[["r.CDRfast_mean"]])

sim_calib <- run_benchmarking(l.inputs1, "", sample1)

sim_calib$prev$fig_prev_by_age
sim_calib$prev$fig_prev_by_raceage

sim_calib$mort$result_plot
sim_calib$reside_time$fig_reside_time
sim_calib$age_onset$fig_age_onset
summary(sim_calib$sim$aggregated_results_totpop$age_at_onset)

# STARTING COHORT ####
tableone_vars <- c("AGE","SEX","RACEETH","EDU","INCOME","MEDBUR","HCARE","APOE4","SYN")
factor_vars <- list(
  SEX  = c("1" = "Male", "2" = "Female"),
  RACEETH   = c("0" = "Non-Hispanic White", "1" = "Non-Hispanic Black", "2" = "Hispanic"),
  INCOME = c("0" = "Low Income", "1" = "Medium Income", "2" = "High Income"),
  HCARE = c("0" = "No regular healthcare provider", "1" = "Has regular healthcare provider"),
  SYN  = c("0" = "Cognitively Healthy", "1" = "Cognitively Impaired"),
  APOE4  = c("0" = "Not Carrier", "1" = "Carrier")
)

cycle1 <- as.data.frame(t(sim_calib$sim$output[1, , ]))  %>% # transpose output so rows are indivs and cols are attrs
  mutate(across(names(factor_vars), ~ factor(.x, levels = names(factor_vars[[cur_column()]]), 
                                             labels = factor_vars[[cur_column()]])))

table1 <- CreateTableOne(vars = tableone_vars, data = cycle1, factorVars = names(factor_vars)) %>%
  print(showAllLevels = TRUE, printToggle = FALSE)

mci_reside_times <- colSums(sim_calib$sim$output[,"SEV",] == 0, na.rm = TRUE)
median_time_mci <- summary(mci_reside_times[mci_reside_times > 0])
dem_reside_times <- colSums(sim_calib$sim$output[,"SEV",] > 0, na.rm = TRUE)
median_time_dem <- summary(dem_reside_times[dem_reside_times > 0])

# Calculate deviation from benchmark prevalence
benchmark_prevalence <- benchmark_prev_by_age %>%
  filter(condition == "dem_manly") %>%
  rename(prev_benchmark = prev)

deviation_prevalence <- sim_calib$prev$dat %>%
  filter(condition == "dem") %>%
  rename(prev_model = prev) %>%
  inner_join(benchmark_prevalence, by = c("age","raceeth")) %>%
  select(age, raceeth, prev_model, prev_benchmark) %>%
  mutate(difference_abs = prev_model - prev_benchmark,
         difference_pct = (prev_model - prev_benchmark) / prev_benchmark)

avg_abs_difference_prevalence <- mean(abs(deviation_prevalence$difference_abs))
avg_abs_difference_prevalence

avg_pct_difference_prevalence <- mean(abs(deviation_prevalence$difference_pct))
avg_pct_difference_prevalence

# Calculate deviation from benchmark mortality
deviation_mortality <- sim_calib$mort$dat %>%
  select(age, model_rate_1000, benchmark_rate_1000, residual, cum_model_rate, cum_benchmark_rate, cum_residual) %>%
  mutate(residual_pct = residual / benchmark_rate_1000,
         cum_residual_pct = cum_residual / cum_benchmark_rate) %>%
  filter(is.finite(residual))

avg_abs_difference_annual_mortality_75under <- mean(abs(deviation_mortality$residual[deviation_mortality$age < 75])) 
avg_abs_difference_annual_mortality_75plus <- mean(abs(deviation_mortality$residual[deviation_mortality$age >= 75])) 

avg_pct_difference_annual_mortality_75under <- mean(abs(deviation_mortality$residual_pct[deviation_mortality$age < 75])) 
avg_pct_difference_annual_mortality_75plus <- mean(abs(deviation_mortality$residual_pct[deviation_mortality$age >= 75])) 

avg_abs_difference_cum_mortality <- mean(abs(deviation_mortality$cum_residual))
avg_pct_difference_cum_mortality <- mean(abs(deviation_mortality$cum_residual_pct))

max_abs_difference_cum_mortality <- deviation_mortality %>%
  filter(cum_residual == max(cum_residual))


# Get dimensions of your simulation output
num_cycles <- dim(sim_calib$sim$output)[1]     # Number of cycles (time steps)
num_individuals <- dim(sim_calib$sim$output)[3] # Number of individuals

# Initialize an empty vector to store the indices of individuals who revert
reverting_individual_indices <- c()

# Loop through each individual
for (i in 1:num_individuals) {
  # Get the SYN status for the current individual across all cycles
  syn_history <- sim_calib$sim$output[, "SYN", i]
  syn_history[is.na(syn_history)] <- 9
  
  # Check for the transition from SYN == 1 in previous cycle to SYN == 0 in current cycle
  # We start from the second cycle (index 2) because we need a 'previous' cycle
  # `syn_history[1:(num_cycles - 1)]` gives previous cycle's SYN
  # `syn_history[2:num_cycles]` gives current cycle's SYN
  has_reverted_in_any_cycle <- any(syn_history[2:num_cycles] == 0 & syn_history[1:(num_cycles - 1)] == 1)
  
  # If this individual has reverted in at least one cycle, add their index to our list
  if (has_reverted_in_any_cycle) {
    reverting_individual_indices <- c(reverting_individual_indices, i)
  }
}

# Print the individual indices
print("Individual indices of those who reverted from SYN 1 to SYN 0:")
if (length(reverting_individual_indices) > 0) {
  print(reverting_individual_indices)
  print(paste("Total number of unique individuals who reverted:", length(reverting_individual_indices)))
} else {
  print("No individuals reverted from SYN 1 to SYN 0 in the simulation.")
}


