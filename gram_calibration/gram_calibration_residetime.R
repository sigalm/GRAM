# ********************************************************************************
# ======= GRAM CALIBRATION ======= 
# ********************************************************************************

# TECHNICAL ##### 
# run only once, do not delete
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
# ***************


# SCENARIOS ####
## Test Run: Natural progression, no treatment ####

# CALIBRATION #####
# Calibrate average rate of change in CDR-SB to achieve expected time in MCI and dementia
# Target durations from Hale 2020 systematic review (https://www.sciencedirect.com/science/article/pii/S2352827319304549)
# Adjusting values for mean rate of change and standard deviation, first for slow, then for fast progressors.

l.inputs_optim <- l.inputs
l.inputs_optim[["r.CDRslow_sd1"]] <- 0
l.inputs_optim[["r.CDRfast_sd1"]] <- 0

# 1. Starting with AD-specific reside times
target_mci_dur <- 3.74
target_dem_dur <- 3.65 + 3.20 + 1.30  # 1.30 is overall severe (all dementias) since an AD-specific value is not available

# MCI 
of_mci <- function(params) {
  print(params)
  l.inputs_optim[["r.CDRslow_mean"]] <- params[1]
  # l.inputs_optim[["r.CDRslow_sd1"]] <- params[2]
  result <- f.wrap_run(l.inputs_optim, printLevel = 0)
  mci_dur <- result$aggregated_results$reside_time$censored[6,2]
  mci_difference <- abs(mci_dur - target_mci_dur)
  return(mci_difference)
}

# Test 1: Using method L-BFGS-B (most common approach for parameters with constraints)
initial_guess_mci <- c(0.3,l.inputs_optim[["r.CDRslow_sd1"]])
result_optim_mci <- optim(par=initial_guess_mci[1], fn = of_mci, method = "L-BFGS-B", lower = 0.1, upper = 2)
result_optim_mci

l1 <- l.inputs
l1[["r.CDRslow_mean"]] <- result_optim_mci$par

t1 <- f.wrap_run(l1)
t1$fig.progression$fig.progression_true
t1$aggregated_results_totpop$MCI.sum
t1$aggregated_results_totpop$reside_time_by_age_of_onset


# Dementia

of_dem <- function(params) {
  print(params)
  l.inputs_optim[["r.CDRfast_mean"]] <- params[1]
  # l.inputs_optim[["r.CDRfast_sd1"]] <- params[2]
  result <- f.wrap_run(l.inputs_optim, printLevel = 0)
  dem_dur <- result$aggregated_results$SEV1.sum + result$aggregated_results$SEV2.sum +result$aggregated_results$SEV3.sum
  dem_difference <- abs(dem_dur - target_dem_dur)
  return(dem_difference)
}

# Test 1: Using method L-BFGS-B (most common approach for parameters with constraints)
initial_guess_dem <- c(l.inputs_optim[["r.CDRfast_mean"]],l.inputs_optim[["r.CDRfast_sd1"]])
result_optim_dem <- optim(par=initial_guess_dem[1], fn = of_mci, method = "L-BFGS-B", lower = 0.5, upper = 2.5)
result_optim_dem

# calculate standard dev around mean duration
l2 <- l.inputs
l2[["r.CDRfast_mean"]] <- result_optim_dem$par

t2 <- f.wrap_run(l1)
t2$fig.progression$fig.progression_true
t2$aggregated_results_totpop$SEV1.sum + t2$aggregated_results_totpop$SEV2.sum + t2$aggregated_results_totpop$SEV3.sum
t2$aggregated_results_totpop$reside_time_by_age_of_onset

t3 <- f.wrap_run(l.inputs)
t3$aggregated_results_totpop$SEV1.sum + t3$aggregated_results_totpop$SEV2.sum + t3$aggregated_results_totpop$SEV3.sum
t3$aggregated_results_totpop$reside_time_by_age_of_onset
