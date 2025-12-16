######################################## GRAM-US PAPER 2: BHA TESTING STRATEGIES ########################################

# Setup, load libraries, source scripts
source("gram_model/gram_setup.r")
source("gram_model/gram_simulation.r")
source("gram_calibration/gram_benchmarking_helpers.r")
source("gram_model/gram_helpers/source_all.r")
library(tableone)
sample1 <- readRDS("gram_data/acs_data/acs_age50_RACE-revised.rds")
# scenario_files <- list.files("gram_model/gram_config/bha_scenarios", full.names = TRUE, pattern = "\\.R$")

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
config_dir <- "gram_analysis/2025_Jul_paper2_testingstrategies/bha_scenarios"
output_dir <- "gram_analysis/2025_Jul_paper2_testingstrategies/sim_results"


## Run scenarios ####
for (scen in names(scenario_list[4])) {
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
  saveRDS(test_data, file = file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
}

# Step 2: generate plots
plot_list <- list()

# all plots have max y axis 100,000 (total N in simulation)
# can change using the y_max argument of plot_test_results()

# # linear y axis - all lines
# for (scen in names(scenario_list[1:3])) {
#   test_data <- readRDS(file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
#   test_plot <- plot_test_results(test_data, ages = 65:80, plot_title = "All Test Results", scenario_name = scen)
#   plot_list[[scen]] <- test_plot
# }
# 
# plot_bhapos_linear <- (plot_list[["u3bhapos"]] + plot_list[["s1bhapos"]] + plot_list[["r1bhapos"]]) + 
#   plot_layout(ncol = 3, guides = "collect") + 
#   plot_annotation(theme = theme(legend.position = "bottom"))
# ggsave("gram_analysis/2025_Jul_paper2_testingstrategies/plots/all-lines-linear.jpeg", plot_bhapos_linear, height = 4, width = 10)  
# 
# 
# # log y axis - all lines
# for (scen in names(scenario_list[1:3])) {
#   test_data <- readRDS(file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
#   test_plot <- plot_test_results(test_data, ages = 65:80, plot_title = "All Test Results", scenario_name = scen, y_transform = "log10")
#   plot_list[[scen]] <- test_plot
# }
# 
# plot_bhapos_log10 <- (plot_list[["u3bhapos"]] + plot_list[["s1bhapos"]] + plot_list[["r1bhapos"]]) + 
#   plot_layout(ncol = 3, guides = "collect") + 
#   plot_annotation(theme = theme(legend.position = "bottom"))
# ggsave("gram_analysis/2025_Jul_paper2_testingstrategies/plots/all-lines-semilog10.jpeg", plot_bhapos_log10, height = 4, width = 10) 
# 

# linear y axis - no early positives (i.e., dotted lines)
for (scen in names(scenario_list[1:3])) {
  test_data <- readRDS(file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
  test_plot <- plot_test_results(test_data, ages = 65:80, plot_title = "All Test Results", scenario_name = scen, show_early_pos = FALSE, y_max = 75000)
  plot_list[[scen]] <- test_plot
}

plot_bhapos_noearlypos <- (plot_list[["u3bhapos"]] + plot_list[["s1bhapos"]] + plot_list[["r1bhapos"]]) + 
  plot_layout(ncol = 3, guides = "collect") + 
  plot_annotation(theme = theme(legend.position = "bottom",
                                legend.box = "vertical", legend.spacing.y = unit(0.1, "cm")))
ggsave("gram_analysis/2025_Jul_paper2_testingstrategies/plots/no-early-positives.jpeg", plot_bhapos_noearlypos, height = 5, width = 8) 


# # linear y axis - no non-testers (i.e., x marked lines)
# for (scen in names(scenario_list[1:3])) {
#   test_data <- readRDS(file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
#   test_plot <- plot_test_results(test_data, ages = 65:80, plot_title = "All Test Results", scenario_name = scen, show_non_testers = FALSE)
#   plot_list[[scen]] <- test_plot
# }
# 
# plot_bhapos_nonontesters <- (plot_list[["u3bhapos"]] + plot_list[["s1bhapos"]] + plot_list[["r1bhapos"]]) + 
#   plot_layout(ncol = 3, guides = "collect") + 
#   plot_annotation(theme = theme(legend.position = "bottom"))
# ggsave("gram_analysis/2025_Jul_paper2_testingstrategies/plots/no-non-testers.jpeg", plot_bhapos_nonontesters, height = 4, width = 10) 
# 
# 
# # linear y axis - no non-testers or early positives
# for (scen in names(scenario_list[1:3])) {
#   test_data <- readRDS(file.path("gram_analysis/2025_Jul_paper2_testingstrategies/test_perf_results", paste0(scen, ".rds")))
#   test_plot <- plot_test_results(test_data, ages = 65:80, plot_title = "All Test Results", scenario_name = scen, show_non_testers = FALSE, show_early_pos = FALSE)
#   plot_list[[scen]] <- test_plot
# }
# 
# plot_bhapos_nonontesters_noearlypositives <- (plot_list[["u3bhapos"]] + plot_list[["s1bhapos"]] + plot_list[["r1bhapos"]]) + 
#   plot_layout(ncol = 3, guides = "collect") + 
#   plot_annotation(theme = theme(legend.position = "bottom"))
# ggsave("gram_analysis/2025_Jul_paper2_testingstrategies/plots/no-non-testers-or-early-positives.jpeg", plot_bhapos_nonontesters_noearlypositives, height = 4, width = 10) 


## Reporting methods ####
# Table 1: Testing likelihood by strategy and cognitive state
reactive_testing_likelihood <- l.inputs_calibrated$m.cogcon_spon[1,-1]
l.inputs_calibrated$m.cogcon_elic
l.inputs_calibrated$m.cogcon

tab1 <- flextable(as.data.frame(rbind(
  l.inputs_calibrated$m.cogcon_spon[1,-1],
  l.inputs_calibrated$m.cogcon_elic[1,-1],
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

