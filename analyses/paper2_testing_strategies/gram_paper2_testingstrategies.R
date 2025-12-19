######################################## GRAM-US PAPER 2: BHA TESTING STRATEGIES ########################################

# Setup, load libraries, source scripts
source("model/setup.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")
source("model/helpers/source_all.R")
library(tableone)
sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")

# Calibration
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 100000)

# List of scenario config files and output names
scenario_list <- list(
  
  u3bhapos = "scenario_u3bhapos_config.R",
  s1bhapos = "scenario_s1bhapos_config.R",
  r1bhapos = "scenario_r1bhapos_config.R",
  
  u3pcppos = "scenario_u3pcppos_config.R",
  s1pcppos = "scenario_s1pcppos_config.R",
  r1pcppos = "scenario_r1pcppos_config.R"
)

# Directory paths
config_dir <- "analyses/paper2_testing_strategies/bha_scenarios"
output_dir <- "analyses/paper2_testing_strategies/sim_results"


## Run scenarios ####
for (scen in names(scenario_list[1:3])) {
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
for (scen in names(scenario_list[1:3])) {
  output <- latest_rds(scen)$output
  test_data <- post_processing_outputs(output)
  saveRDS(test_data, file = file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, ".rds")))
}

# Step 2: generate plots
test_data_combined <- data.frame()
for (i in seq_along(names(scenario_list[1:3]))) {
  scen <- names(scenario_list)[i]
  test_data <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, ".rds"))) %>%
    mutate(scenario = scen)
  test_data_combined <- rbind(test_data_combined, test_data)
}

subtitles <- c("Inclusive testing, every 3 years",
               "Selective testing, annual",
               "Reactive testing, annual")
names(subtitles) <- names(scenario_list)[1:3]
plot_test_results(test_data_combined, ages = 65:80, show_early_pos = FALSE, scenario_names = subtitles, y_max = 75000)

ggsave("analyses/paper2_testing_strategies/plots/no-early-positives.jpeg", height = 10, width = 8) 



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
u3bhapos <- latest_rds("u3bhapos")$output
subset_age_65 <- as.data.frame(t(u3bhapos[65-50+1,,])) # rows are IDs, cols are attributes

prev_in_undx_65 <- subset_age_65 %>%
  filter(ALIVE == 1, DX == 0) %>%
  mutate(status = case_when(SYN < 1 ~ "h",
                            SEV == 0 ~ "mci",
                            SEV >= 1 ~ "dem")) %>%
  group_by(status) %>%
  summarise(n = n()) %>%
  mutate(prev = n/sum(n))















## Early diagnoses ####

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

