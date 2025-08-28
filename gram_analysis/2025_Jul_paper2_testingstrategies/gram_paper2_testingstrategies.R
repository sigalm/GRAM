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
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 100000)

# List of scenario config files and output names
scenario_list <- list(
  u1      = "scenario_u1_config.R",
  s1      = "scenario_s1_config.R",
  r1      = "scenario_r1_config.R",
  #  u2      = "scenario_u2_config.R",
  u1bhapos = "scenario_u1bhapos_config.R",
  s1bhapos = "scenario_s1bhapos_config.R",
  r1bhapos = "scenario_r1bhapos_config.R",
  #  u1pcppos = "scenario_u1pcppos_config.R",
  u1_bhags = "scenario_u1_bhags_config.R",
  s1_bhags = "scenario_s1_bhags_config.R",
  r1_bhags = "scenario_r1_bhags_config.R"
)

# Directory paths
config_dir <- "gram_model/gram_config/bha_scenarios"
output_dir <- "gram_analysis/2025_Jul_paper2_testingstrategies/sim_results"

for (scen in names(scenario_list[8:9])) {
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


# Load results
u1 <- readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_u1_sim.rds")
s1 <- readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_s1_sim.rds")
r1 <- readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_r1_sim.rds")

u1bhapos <- readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_u1bhapos_sim_20250827_152049.rds")
s1bhapos <-  readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_s1bhapos_sim_20250827_152734.rds")
r1bhapos <-  readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_r1bhapos_sim_20250827_153350.rds")

u1_bhags <- readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_u1_bhags_sim_20250827_143819.rds")
s1_bhags <-  readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_s1_bhags_sim_20250827_162918.rds")
r1_bhags <-  readRDS("gram_analysis/2025_Jul_paper2_testingstrategies/sim_results/scenario_r1_bhags_sim_20250827_163405.rds")

View(plot_results(
  list(
    #     config_r1, config_p1, 
    config_u1bhapos, config_u1pcppos), 
  #   list(
  #     # r1, p1, 
  u1bhapos, u1pcppos))$results_table
# )

scenario_colors <- c(
  R1 = "aquamarine3", R1BHAPOS = "aquamarine3", R1_BHAGS = "aquamarine3",
  S1 = "darksalmon", S1BHAPOS = "darksalmon", S1_BHAGS = "darksalmon",
  U1 = "darkgrey", U2 = "darkgrey", U1BHAPOS = "darkgrey", U1PCPPOS = "darkgrey", U1_BHAGS = "darkgrey"
)
scenario_linetypes <- c(
  R1 = "solid", R1BHAPOS = "dashed", R1_BHAGS = "dotted",
  S1 = "solid", S1BHAPOS = "dashed", S1_BHAGS = "dotted",
  U1 = "solid", U2 = "solid", U1BHAPOS = "dashed", U1PCPPOS = "dotted",
  U1_BHAGS = "dotdash"
) 

######## Pie charts ####
make_gram_pie_chart <- function(output_array, cycles = 21:30, plot_title = "Cumulative Testing Status") {
  # Extract BHA and SYN matrices for the specified cycles
  bha_mat <- output_array[cycles, "BHA", ]
  syn_mat <- output_array[cycles, "SYN", ]
  
  # Determine eligibility and testing status
  is_eligible <- apply(bha_mat, 2, function(x) any(!is.na(x) & (x != -9)))
  ever_tested <- apply(bha_mat, 2, function(x) any(x %in% c(0, 1), na.rm = TRUE))
  never_tested <- is_eligible & !ever_tested
  
  never_tested_idx <- which(never_tested)
  never_tested_healthy <- 0
  never_tested_impaired <- 0
  
  for (i in never_tested_idx) {
    bha_i <- bha_mat[, i]
    syn_i <- syn_mat[, i]
    eligible_cycles <- which(!is.na(bha_i) & (bha_i != -9))
    syn_eligible <- syn_i[eligible_cycles]
    
    if(length(syn_eligible) == 0) next
    if (any(syn_eligible == 1, na.rm = TRUE)) {
      never_tested_impaired <- never_tested_impaired + 1     # If ever impaired, count as impaired (false negative)
    } else if (all(syn_eligible %in% c(0, 0.5), na.rm = TRUE)) {
      never_tested_healthy <- never_tested_healthy + 1       # If never impaired, count as healthy (true negative)
    }
  }
  
  tested_idx <- which(is_eligible & ever_tested)
  TP <- 0
  FP <- 0
  TN <- 0
  FN <- 0
  
  consistently_tp_tn <- 0
  fp_to_tn <- 0
  fn_to_tp <- 0
  
  for (i in tested_idx) {
    bha_i <- bha_mat[, i]
    syn_i <- syn_mat[, i]
    eligible_cycles <- which(!is.na(bha_i) & (bha_i != -9))
    bha_eligible <- bha_i[eligible_cycles]
    syn_eligible <- syn_i[eligible_cycles]
    # Only consider cycles where a test was actually considered (BHA==0 or 1)
    tested_cycles <- which(bha_eligible %in% c(0, 1))
    if (length(tested_cycles) == 0) next
    bha_tested <- bha_eligible[tested_cycles]
    syn_tested <- syn_eligible[tested_cycles]
    ever_pos <- any(bha_tested == 1, na.rm = TRUE)
    ever_neg <- any(bha_tested == 0, na.rm = TRUE)
    ever_imp <- any(syn_tested == 1, na.rm = TRUE)
    always_healthy <- all(syn_tested %in% c(0, 0.5), na.rm = TRUE)
    
    if (ever_pos && ever_imp) {
      TP <- TP + 1
    } else if (ever_pos && always_healthy) {
      FP <- FP + 1
    } else if (ever_neg && !ever_pos && always_healthy) {
      TN <- TN + 1
    } else if (ever_neg && !ever_pos && ever_imp) {
      FN <- FN + 1
    }
    # If both pos and neg and both healthy and impaired, will be counted as TP (priority: TP > FP > FN > TN)
    # Does not require contemporaneous identification. Rather, counts if individuals were ever "caught" by the system

    # Test trajectory
    status_vec <- rep(NA_character_, length(bha_tested))
    status_vec[bha_tested == 1 & syn_tested == 1] <- "TP"
    status_vec[bha_tested == 1 & syn_tested %in% c(0, 0.5)] <- "FP"
    status_vec[bha_tested == 0 & syn_tested %in% c(0, 0.5)] <- "TN"
    status_vec[bha_tested == 0 & syn_tested == 1] <- "FN"
    status_vec <- status_vec[!is.na(status_vec)]
    
    # Consistently TP or TN
    if (length(status_vec) > 0 && all(status_vec %in% c("TP", "TN"))) {
      consistently_tp_tn <- consistently_tp_tn + 1
    }
    # FP corrected to TN: at least one FP, and a TN occurs after the first FP
    fp_idx <- which(status_vec == "FP")
    tn_idx <- which(status_vec == "TN")
    if (length(fp_idx) > 0 && length(tn_idx) > 0 && any(tn_idx > min(fp_idx))) {
      fp_to_tn <- fp_to_tn + 1
    }
    # FN corrected to TP: at least one FN, and a TP occurs after the first FN
    fn_idx <- which(status_vec == "FN")
    tp_idx <- which(status_vec == "TP")
    if (length(fn_idx) > 0 && length(tp_idx) > 0 && any(tp_idx > min(fn_idx))) {
      fn_to_tp <- fn_to_tp + 1
    }
  }
  
  counts <- c(
    "True Positive" = TP,
    "False Positive" = FP,
    "True Negative" = TN,
    "False Negative" = FN,
    "Never tested & Healthy" = never_tested_healthy,
    "Never tested & Impaired" = never_tested_impaired
  )
  
  bha_df_pies <- data.frame(status = names(counts), count = as.numeric(counts))
  
  piechart <- ggplot(bha_df_pies, aes(x = "", y = count, fill = status)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(title = plot_title, fill = "") +
    theme_void() +
    scale_fill_manual(values = c("red", "pink", "skyblue", "orange", "forestgreen", "goldenrod"))
  
  test_trajectory <- data.frame(
    consistently_tp_tn = consistently_tp_tn,
    fp_to_tn = fp_to_tn,
    fn_to_tp = fn_to_tp
  )
  return(list(df = bha_df_pies, piechart = piechart, test_trajectory_stats = test_trajectory))
}

u1_pie <- make_gram_pie_chart(output_array = u1$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Universal Testing")
u1_pie

s1_pie <- make_gram_pie_chart(output_array = s1$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Selective Testing")
s1_pie

r1_pie <- make_gram_pie_chart(output_array = r1$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Reactive Testing")
r1_pie

u1bhapos_pie <- make_gram_pie_chart(output_array = u1bhapos$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Universal Testing Until BHA+")
u1bhapos_pie

s1bhapos_pie <- make_gram_pie_chart(output_array = s1bhapos$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Selective Testing Until BHA+")
s1bhapos_pie

r1bhapos_pie <- make_gram_pie_chart(output_array = r1bhapos$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Reactive Testing Until BHA+")
r1bhapos_pie

u1_bhags_pie <- make_gram_pie_chart(output_array = u1_bhags$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Universal Testing with Full BHA")
u1_bhags_pie

s1_bhags_pie <- make_gram_pie_chart(output_array = s1_bhags$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Selective Testing with Full BHA")
s1_bhags_pie

r1_bhags_pie <- make_gram_pie_chart(output_array = r1_bhags$output, cycles = 21:30, plot_title = "Cumulative Testing Status, Ages 70-79, Reactive Testing with Full BHA")
r1_bhags_pie



#########
plot1_nostop <- plot_results(
  list(
    r1,
    s1,
    u1, u1_bhags))

plot1_nostop$results_figure

plot2_universals <- plot_results(
  list(
    u1, u2,
    u1bhapos, u1pcppos))$results_figure

x <- plot_results(
  list(
    r1, r1bhapos,
    s1, s1bhapos,
    u1, u2,
    u1bhapos, u1pcppos))$results_figure


plot_results_cumulative(
  list(
    r1,
    p1,
    u1))$results_figure

plot_results_cumulative(
  list(
    config_u1, config_u2,
    config_u1bhapos, config_u1pcppos), 
  list(
    u1, u2,
    u1bhapos, u1pcppos))$results_figure

plot_results_cumulative(
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
