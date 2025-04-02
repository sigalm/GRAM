# ********************************************************************************
# ======= GRAM PAPER 1: BHA TEST PERFORMANCE AND FALSE POSITIVE RATE ======= 
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
l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(1, 1, 1) 
l.inputs1[["r.CDRslow_mean"]] <-  (seq(0.25, 0.75, length.out = 51)^2) * (2 * l.inputs[["r.CDRslow_mean"]])
l.inputs1[["p.HCARE_start"]] <- c(0.25,0.75)

sim_calib <- run_benchmarking(l.inputs1, "calibration run", sample1)

# STARTING COHORT ####
tableone_vars <- c("AGE","SEX","RACEETH","EDU","INCOME","MEDBUR","APOE4","SYN")
factor_vars <- list(
  SEX  = c("1" = "Male", "2" = "Female"),
  RACEETH   = c("0" = "Non-Hispanic White", "1" = "Non-Hispanic Black", "2" = "Hispanic"),
  SYN  = c("0" = "Cognitively Healthy", "1" = "Cognitively Impaired"),
  APOE4  = c("0" = "Not Carrier", "1" = "Carrier"),
  INCOME = c("0" = "Low Income", "1" = "Medium Income", "2" = "High Income")
)

cycle1 <- as.data.frame(t(sim_calib$sim$output[1, , ]))  %>% # transpose output so rows are indivs and cols are attrs
  mutate(across(names(factor_vars), ~ factor(.x, levels = names(factor_vars[[cur_column()]]), 
                                             labels = factor_vars[[cur_column()]])))

table1 <- CreateTableOne(vars = tableone_vars, data = cycle1, factorVars = names(factor_vars)) %>%
  print(showAllLevels = TRUE, printToggle = FALSE)



# SCENARIOS ####

## BHA at age 70 ####
###  Scenario 1a : Screen everyone with regular healthcare provider (HCARE) without prior diagnosis ####
l.inputs1a <- l.inputs1
l.inputs1a[["scenario"]] <- "Scenario 1a: Universal BHA" 
# l.inputs1a[["m.cogcon"]]  # All probs are 1, so everyone gets BHA regardless of cognitive concerns status
l.inputs1a[["sens_BHA"]] 
l.inputs1a[["spec_BHA"]] 

scenario1a <- f.wrap_run(l.inputs1a, microdata = sample1)

### Scenario 1b: Screen those with regular healthcare provider (HCARE) without prior diagnosis AND systematic query of cognitive concerns ####
l.inputs1b <- l.inputs1
l.inputs1b[["scenario"]] <- "Scenario 1b: Systematic query in AWV or equivalent"
l.inputs1b[["m.cogcon"]] <- l.inputs[["m.cogcon_elic"]]
l.inputs1b[["sens_BHA"]] <- c(0.80, 1.00)
l.inputs1b[["spec_BHA"]] <- 0.72

scenario1b <- f.wrap_run(l.inputs1b, microdata = sample1)

### Scenario 1c: Screen those with regular healthcare provider (HCARE) without prior diagnosis AND spontaneously endorses cognitive concerns ####
l.inputs1c <- l.inputs1
l.inputs1c[["scenario"]] <- "Scenario 1c: Spontaneous endorsement of cognitive concerns"
l.inputs1c[["m.cogcon"]] <- l.inputs[["m.cogcon_spon"]]
l.inputs1c[["sens_BHA"]] <- c(0.87, 1.00)
l.inputs1c[["spec_BHA"]] <- 0.60

scenario1c <- f.wrap_run(l.inputs1c, microdata = sample1)



# Present results

compute_results <- function(inputs, scenario) {
  map_dfr(50:100, function(age) {
    testperf <- f.analyze_test_performance(scenario, age - 50 + 1)
    
    tibble(
      Scenario = inputs[["scenario"]],
      Age = age,
      N_Tests = testperf$total_tests,
      TP = sum(testperf$concordance_table[c("mci", "dem"), "BHApos"]),
      FN = sum(testperf$concordance_table[c("mci", "dem"), "BHAneg"]),
      TN = testperf$concordance_table["h", "BHAneg"],
      FP = testperf$concordance_table["h", "BHApos"],
      PPV = round(testperf$ppv, 3),
      NPV = round(testperf$npv, 3)
    )
  })
}

inputs_list <- list(l.inputs1a, l.inputs1b, l.inputs1c)
scenario_list <- list(scenario1a, scenario1b, scenario1c)

all_results <- map2_dfr(inputs_list, scenario_list, compute_results)
# flextable(all_results)

# results_70 <- all_results %>%
#   filter(Age == 70)

# Get prevalence of impairment in those without prior diagnosis at age 70 ####

get_prev_no_dx <- function(l.inputs, output) {
  n_cycles <- l.inputs[["n.cycle"]]
  prev_no_dx <- matrix(NA, nrow = n_cycles, ncol = 5, 
                       dimnames = list(NULL, c("healthy", "mci", "mil", "mod", "sev")))
  
  # Identify individuals undiagnosed in the previous cycle (cycle 2 onward)
  prior_no_dx <- rbind(NA, output[-n_cycles, "DX", ] == 0)  # Shift DX to align with previous cycle
  alive_now <- output[, "ALIVE", ] == 1  # Alive in the current cycle
  
  # Number of people alive without prior diagnosis
  n_no_dx <- rowSums(prior_no_dx & alive_now, na.rm = TRUE)
  
  # Compute prevalence for each category
  prev_no_dx[, "healthy"] <- rowSums(output[, "SYN", ] == 0 & prior_no_dx & alive_now, na.rm = TRUE) / n_no_dx
  prev_no_dx[, "mci"] <- rowSums(output[, "SEV", ] == 0 & prior_no_dx & alive_now, na.rm = TRUE) / n_no_dx
  prev_no_dx[, "mil"] <- rowSums(output[, "SEV", ] == 1 & prior_no_dx & alive_now, na.rm = TRUE) / n_no_dx
  prev_no_dx[, "mod"] <- rowSums(output[, "SEV", ] == 2 & prior_no_dx & alive_now, na.rm = TRUE) / n_no_dx
  prev_no_dx[, "sev"] <- rowSums(output[, "SEV", ] == 3 & prior_no_dx & alive_now, na.rm = TRUE) / n_no_dx
  
  return(as.data.frame(prev_no_dx))
}
prev_no_dx_1a <- get_prev_no_dx(l.inputs1a, scenario1a$output)
prev_no_dx_1a[21,"mci"]
sum(prev_no_dx_1a[21,c("mil","mod","sev")])



# FIGURES ####
pv_results <- all_results %>%
  pivot_longer(cols = c(PPV, NPV), names_to = "Metric", values_to = "Value") %>%
  select(Scenario, Age, Metric, Value) %>%
  filter(Age <= 90)   # Older ages have small sample problems, so remove

fig_results <- ggplot(data = pv_results, aes(x = Age, y = Value, color = Scenario)) +
  geom_line() + 
  facet_wrap(~ Metric) + 
  theme_minimal()
  
fig_results




# Export figures (change file path extension if you want a different file type e.g., pdf, jpeg)
# ggsave("[replace with file path].png", plot = fig_scenario1, width = 10, height = 6, dpi = 300)

## Export tables (you can add multiple tables to a single word doc by including multiple tables in the function call
# save_as_docx("[replace with table title]" = reside_time_scenario1, path = "[replace with file path]")

