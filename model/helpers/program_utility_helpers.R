######################################## GRAM HELPER FUNCTIONS: PROGRAM UTILITY HELPERS ########################################
# This script defines program utility helpers.

# Load configuration files
load_scenario <- function(config_file, base_inputs) {
  env <- new.env()
  source(config_file, local = env)
  base_inputs[["scenario"]] <- modifyList(base_inputs[["scenario"]], env$scenario_inputs)
  base_inputs
}


# Read the results of last saved run
latest_rds <- function(prefix, dir = output_dir) {
  files <- list.files(output_dir, pattern = paste0("scenario_", prefix, "_sim_.*\\.rds$"), full.names = TRUE)
  
  if (length(files) == 0) stop("No files found for prefix: ", prefix)
  latest <- sort(files, decreasing = TRUE)[1]
  
  message("Loading: ", latest)
  readRDS(latest)
}