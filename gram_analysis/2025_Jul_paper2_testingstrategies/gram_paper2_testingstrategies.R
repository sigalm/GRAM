######################################## GRAM-US PAPER 2: BHA TESTING STRATEGIES ########################################

# Setup, load libraries, source scripts
source("gram_model/gram_setup.r")
source("gram_model/gram_simulation.r")
source("gram_calibration/gram_benchmarking_helpers.r")
source("gram_model/gram_helpers/source_all.r")
library(tableone)
sample1 <- readRDS("gram_data/acs_data/acs_age50.rds")
# scenario_files <- list.files("gram_model/gram_config/bha_scenarios", full.names = TRUE, pattern = "\\.R$")

# Calibration
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 10000)

config_u1 <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_u1_config.R", l.inputs_calibrated)
config_p1 <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_p1_config.R", l.inputs_calibrated)
config_r1 <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_r1_config.R", l.inputs_calibrated)
config_u2 <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_u2_config.R", l.inputs_calibrated)
config_u1bhapos <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_u1bhapos_config.R", l.inputs_calibrated)
config_p1bhapos <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_p1bhapos_config.R", l.inputs_calibrated)
config_r1bhapos <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_r1bhapos_config.R", l.inputs_calibrated)
config_u1pcppos <- load_scenario("gram_model/gram_config/bha_scenarios/scenario_u1pcppos_config.R", l.inputs_calibrated)


config_u1[["scenario"]][["title"]]     # Check to see right scenario
u1 <- f.wrap_run(config_u1, microdata = sample1)

config_p1[["scenario"]][["title"]]
p1 <- f.wrap_run(config_p1, microdata = sample1)

config_r1[["scenario"]][["title"]]
r1 <- f.wrap_run(config_r1, microdata = sample1)

config_u2[["scenario"]][["title"]]    
u2 <- f.wrap_run(config_u2, microdata = sample1)

config_u1bhapos[["scenario"]][["title"]]    
u1bhapos <- f.wrap_run(config_u1bhapos, microdata = sample1)

config_p1bhapos[["scenario"]][["title"]]    
p1bhapos <- f.wrap_run(config_p1bhapos, microdata = sample1)

config_r1bhapos[["scenario"]][["title"]]    
r1bhapos <- f.wrap_run(config_r1bhapos, microdata = sample1)

config_u1pcppos[["scenario"]][["title"]]
u1pcppos <- f.wrap_run(config_u1pcppos, microdata = sample1)


 View(plot_results(
   list(
#     config_r1, config_p1, 
     config_u1bhapos, config_u1pcppos), 
#   list(
#     # r1, p1, 
     u1bhapos, u1pcppos))$results_table
# )

plot_results(
  list(
    config_r1, config_r1bhapos,
    config_p1, config_p1bhapos,
    config_u1, config_u2,
    config_u1bhapos, config_u1pcppos), 
  list(
    r1, r1bhapos,
    p1, p1bhapos,
    u1, u2,
    u1bhapos, u1pcppos))$results_figure


