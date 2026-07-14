######################################## GRAM-US PAPER 2: BHA TESTING STRATEGIES ########################################

# Setup, load libraries, source scripts
source("model/setup.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")
source("model/helpers/source_all.R")
source("analyses/paper2_testing_strategies/test_performance_helpers.R")
library(tableone)
library(flextable)
sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")

# Calibration
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 100000)

# Analysis plan:
# Strategies analyzed will follow the GRAM-ish one-time testing paper
#   (1) reactive testing
#   (2) selective with eRADAR/EHR-based
#   (3) inclusive 50% random opt-in

# Outcomes reported:
#   (a) overall sensitivity/specificity
#   (b) PPV and NPV 
#   all at two time points: at first test, and after 10 years

# Sensitivity analyses
#   (1) a question-based selective option
#   (2) inclusive with non-random selection (NB, unlike GRAM-ish, total testing not fixed at 50% here)
#   (3) threshold inclusive with random selection?? TBD


# List of scenario config files and output names
scenario_list <- list(
  
  u3bhapos_rand50 = "scenario_u3bhapos_rand50_config.R",
  s1bhapos_emr = "scenario_s1bhapos_emr_config.R",
  r1bhapos = "scenario_r1bhapos_config.R",
  
  s1bhapos_question = "scenario_s1bhapos_question_config.R",
  u3bhapos_nonrand = "scenario_u3bhapos_nonrand_config.R"
  )



# Directory paths
config_dir <- "analyses/paper2_testing_strategies/bha_scenarios"
output_dir <- "analyses/paper2_testing_strategies/sim_results"


## Run scenarios ####
for (scen in names(scenario_list[1:5])) {
  local({
    config_file <- file.path(config_dir, scenario_list[[scen]])
    cat("Running scenario:", scen, "\n")
    config <- load_scenario(config_file, l.inputs_calibrated)
    result <- f.wrap_run(config, microdata = sample1)
    datetime_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
    saveRDS(result, file = file.path(output_dir, paste0("scenario_", scen, "_sim_", datetime_suffix, ".rds")))
    
    rm(config, result)
    invisible(gc())
  })
}


######## Evolution Charts ####
# Step 1: create and save test performance data
for (scen in names(scenario_list[1:5])) {
  output <- latest_rds(scen)$output
  test_data <- post_processing_outputs(output)
  saveRDS(test_data, file = file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, "230126.rds")))
}

# Step 2: generate plots
main_test_data <- data.frame()
for (i in seq_along(names(scenario_list[1:3]))) {
  scen <- names(scenario_list)[i]
  temp_test_data <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, "230126.rds"))) %>%
    mutate(scenario = scen)
  main_test_data <- rbind(main_test_data, temp_test_data)
}

subtitles <- c("Inclusive",
               "Selective",
               "Reactive")
names(subtitles) <- names(scenario_list)[1:3]
p1 <- plot_test_results(main_test_data, ages = 65:80, show_early_pos = FALSE, scenario_names = subtitles, y_max = 75000)
p1



p1 <- p1 + theme(
  text         = element_text(size = 18),      # base size for all text
  axis.title   = element_text(size = 20, face = "bold"),
  axis.text    = element_text(size = 16),
  legend.text  = element_text(size = 16),
  legend.title = element_text(size = 18, face = "bold"),
  strip.text   = element_text(size = 20, face = "bold"),  # "Reactive/Selective/Inclusive" labels
  plot.title   = element_text(size = 24, face = "bold"),
  legend.key.size = unit(1.2, "cm")
)

ggsave("analyses/paper2_testing_strategies/plots/no-early-positives_reordered.jpeg",
       plot = p1,
       height = 10, width = 14,   # wider to give 3 panels more breathing room
       dpi = 300)

plot_testers(main_test_data, 
             ages = 65:80, 
             scenario_names = subtitles) 

ggsave("analyses/paper2_testing_strategies/plots/no-early-positives_reordered.jpeg", plot = p1, height = 10, width = 8) 

## Sensitivity analyses ####
### Question-based selective strategy
scen <- names(scenario_list)[4]
test_data_question_selective <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, "230126.rds"))) %>%
  mutate(scenario = scen)
test_data_selective <- rbind(main_test_data, test_data_question_selective) %>%
  filter(scenario %in% names(scenario_list[c(2,4)]))
plot_test_results(test_data_selective, ages = 65:80, show_early_pos = FALSE, y_max = 75000)

### Non-random inclusive
scen <- names(scenario_list)[5]
test_data_nonrand_inclusive <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, "230126.rds"))) %>%
  mutate(scenario = scen)
test_data_inclusive <- rbind(main_test_data, test_data_nonrand_inclusive) %>%
  filter(scenario %in% names(scenario_list[c(1,5)]))
plot_test_results(test_data_inclusive, ages = 65:80, show_early_pos = FALSE, y_max = 75000)



## Reporting methods ####
# Table 1: Testing likelihood by strategy and cognitive state
reactive_testing_likelihood <- l.inputs_calibrated$m.cogcon_reactive[1,-1]
l.inputs_calibrated$m.cogcon_selective
l.inputs_calibrated$m.cogcon

tab1 <- flextable(as.data.frame(rbind(
  l.inputs_calibrated$m.cogcon_reactive[1,-1],
  l.inputs_calibrated$m.cogcon_selective[1,-1],
  l.inputs_calibrated$m.cogcon[1,-1]
)))

l.inputs_calibrated$sens_BHAGS
l.inputs_calibrated$spec_BHAGS
l.inputs_calibrated$rr.cogcon_prior

## Reporting results ####
u3bhapos <- latest_rds("u3bhapos_rand50")$output
s1bhapos <- latest_rds("s1bhapos_emr")$output
r1bhapos <- latest_rds("r1bhapos")$output

subset_age_65 <- as.data.frame(t(u3bhapos[65-50+1,,])) # rows are IDs, cols are attributes
subset_age_67 <- as.data.frame(t(u3bhapos[67-50+1,,])) # rows are IDs, cols are attributes

prev_in_undx_65 <- subset_age_65 %>%
  filter(ALIVE == 1, DX == 0) %>%
  mutate(status = case_when(SYN < 1 ~ "h",
                            SEV == 0 ~ "mci",
                            SEV >= 1 ~ "dem")) %>%
  group_by(status) %>%
  summarise(n = n()) %>%
  mutate(prev = n/sum(n))


prev_in_undx_67 <- subset_age_67 %>%
  filter(ALIVE == 1, DX == 0) %>%
  mutate(status = case_when(SYN < 1 ~ "h",
                            SEV == 0 ~ "mci",
                            SEV >= 1 ~ "dem")) %>%
  group_by(status) %>%
  summarise(n = n()) %>%
  mutate(prev = n/sum(n))


# List of scenario arrays
scenario_arrays <- list(
  u3bhapos = u3bhapos,
  s1bhapos = s1bhapos,
  r1bhapos = r1bhapos
)

# Function to compute prevalence at first test
get_prev_at_first_test <- function(scenario_array) {
  bha_matrix <- scenario_array[, "BHA", ]
  first_bha_cycle <- apply(bha_matrix, 2, function(x) which(x == 0 | x == 1 )[1])
  
  subset_first_test <- t(sapply(seq_along(first_bha_cycle), function(i) {
    cycle <- first_bha_cycle[i]
    if (!is.na(cycle)) {
      scenario_array[cycle, , i]
    } else {
      rep(NA, dim(scenario_array)[2])
    }
  }))
  colnames(subset_first_test) <- dimnames(scenario_array)[[2]]
  subset_first_test <- as.data.frame(subset_first_test)
  
  prev_by_status <- subset_first_test %>%
    filter(ALIVE == 1, DX == 0) %>%
    mutate(status = case_when(SYN < 1 ~ "h",
                              SEV == 0 ~ "mci",
                              SEV >= 1 ~ "dem")) %>%
    group_by(status) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(prev = n / sum(n))
  
  n_tested <- sum(!is.na(first_bha_cycle))
  n_tested_alive <- sum(subset_first_test$ALIVE == 1, na.rm = TRUE)
  prop_alive_at_test <- n_tested_alive / n_tested
  
  list(
    prev_by_status = prev_by_status,
    n_tested = n_tested
  )
}

# Apply to all scenarios
prev_in_first_test <- lapply(scenario_arrays, get_prev_at_first_test)

# Calculate % of undx tested
prev_in_first_test$u3bhapos$n_tested / sum(prev_in_undx_67$n)
prev_in_first_test$s1bhapos$n_tested / sum(prev_in_undx_65$n)
prev_in_first_test$r1bhapos$n_tested / sum(prev_in_undx_65$n)

# Calculate total CI prev at first test
sum(prev_in_first_test$u3bhapos$prev_by_status$prev[c(1,3)])
sum(prev_in_first_test$s1bhapos$prev_by_status$prev[c(1,3)])
sum(prev_in_first_test$r1bhapos$prev_by_status$prev[c(1,3)])


# Make table with result numbers
results_table <- main_test_data %>%
  filter(age %in% c(65,67,75)) %>%
  mutate(age = age,
         scenario = scenario,
         tp = tp + converted_tp,
         fp = fp + early_pos,
         tn = tn + notest_tn,
         fn = fn + notest_fn,
         dead = death,
         .keep = "none") %>%
  flextable()




predictive_value <- main_test_data %>%
  filter(age %in% c(65, 67, 70, 75, 80)) %>%
  mutate(ppv = (tp + converted_tp) / (tp + early_pos + converted_tp + fp),
         npv = ((tn + notest_tn) / (tn + fn + notest_tn + notest_fn)),
         sens = (tp + converted_tp) / (tp + converted_tp + fn + notest_fn),
         spec = (tn + notest_tn) / (tn + notest_tn + early_pos + fp),
         acc = (tp + converted_tp + tn + notest_tn) / (tp + converted_tp + tn + notest_tn + fp + early_pos + fn + notest_fn))

predictive_value_emr <- test_data_emr_selective %>%
  filter(age %in% c(65, 70, 75, 80)) %>%
  mutate(ppv = (tp + converted_tp) / (tp + early_pos + converted_tp + fp),
         npv = ((tn + notest_tn) / (tn + fn + notest_tn + notest_fn)),
         sens = (tp + converted_tp) / (tp + converted_tp + fn + notest_fn),
         spec = (tn + notest_tn) / (tn + notest_tn + early_pos + fp),
         acc = (tp + converted_tp + tn + notest_tn) / (tp + converted_tp + tn + notest_tn + fp + early_pos + fn + notest_fn))



results_table <- data.frame(
  Time = c(rep("First Test",3), rep("By Year 10", 3)),
  Strategy = rep(c("Inclusive", "Selective", "Reactive"), 2),
  Sensitivity = round(predictive_value$sens[c(2,6,11, 4,9,14)] * 100, digits = 1),
  Specificity = round(predictive_value$spec[c(2,6,11, 4,9,14)] * 100, digits = 1),
  PPV = round(predictive_value$ppv[c(2,6,11, 4,9,14)] * 100, digits = 1),
  NPV = round(predictive_value$npv[c(2,6,11, 4,9,14)] * 100, digits = 1),
  Accuracy = round(predictive_value$acc[c(2,6,11, 4,9,14)] * 100, digits = 1)
)


flextable(results_table)

## OLD CODE -- Early diagnoses ####

u1pcppos$output <- add_early_dx(u1pcppos$output, dx_var = "BHA-DX")


# identified at MCI vs at dementia
# hybrid follow up scenarios 
# start writing out key assumptions 

# could do a version with alives only wihtout the red line
p_cum_earlyID_u1pcppos <- plot_cumulative_count(u1pcppos$output, 
                                                variables = list("SYN" = list(variable_name = "SYN", condition_value = 1),
                                                                 "DX" = list(variable_name = "DX", condition_value = 1),
                                                                 "early_dx_1" = list(variable_name = "early_dx_1", condition_value = 1),
                                                                 "early_dx_2" = list(variable_name = "early_dx_2", condition_value = 1),
                                                                 "early_dx_3" = list(variable_name = "early_dx_3", condition_value = 1),
                                                                 "early_dx_4" = list(variable_name = "early_dx_4", condition_value = 1),
                                                                 "early_dx_5" = list(variable_name = "early_dx_5", condition_value = 1),
                                                ),
                                                plot_title = "Cumulative Early Diagnoses",
                                                scenario_name = "u1pcppos")



(p_cum_earlyID_u1pcppos <- plot_cumulative_count(u1pcppos$output, 
                                                 variables = list("SYN" = list(variable_name = "SYN", condition_value = 1))))











# u1bhapos_evolution_plot <- plot_simulation_evolution(u1bhapos$output, plot_title = "Evolution for u1bhapos Scenario")
# print(u1bhapos_evolution_plot)

u1bhapos_evolution_lines_plot <- plot_simulation_evolution_lines(u1bhapos$output, plot_title = "Evolution for u1bhapos Scenario (Lines)")
print(u1bhapos_evolution_lines_plot)

