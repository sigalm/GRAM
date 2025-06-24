# ============== CDR test performance assessment and calibration ========= #

test_performance <- function(TP, FP, TN, FN) {
  
  epsilon <- 1e-10
  
  sens <- TP / (TP + FN + epsilon)
  spec <- TN / (TN + FP + epsilon)
  ppv <- TP / (TP + FP + epsilon)
  npv <- TN / (TN + FN + epsilon)
  
  average <- c("Sens" = mean(sens), "Spec" = mean(spec), 
           "PPV" = mean(ppv), "NPV" = mean(npv))
  
  return(list(
    sens = sens,
    spec = spec,
    ppv = ppv,
    npv = npv,
    average = average
  ))
  
}


# True states
# True severity is coded under "SEV", 0 = MCI, 1 = Mild dementia, 2 = Moderate dementia, 3 = Severe dementia
# In "SEV", both healthy and dead are coded as NA
# Find healthy people in "SYN" (which is true healthy vs. impairment), assign them -9 in "SEV"

true_state <- model$output[-1,"SEV",]
thealthy <- model$output[-1,"SYN",] == 0
true_state[thealthy] <- -9

# Observed states
# Observed severity is coded under "SEV_obs", 0 = MCI, 1 = Mild dementia, 2 = Moderate dementia, 3 = Severe dementia
# In "SEV_obs", both healthy and dead are coded as NA
# Find healthy people in "BHA" (which is observed healthy vs. impairment), assign them -9 in "SEV_obs"

ohealthy <- model$output[-1,"BHA",] == 0
obs_state <- model$output[-1,"SEV_obs",]
obs_state[ohealthy] <- -9


# CDR test performance for ANY IMPAIRMENT vs HEALTHY
# (positive means any impairment i.e., state NOT EQUAL TO -9)
TP_any <- rowSums((obs_state != -9) & (true_state != -9), na.rm = TRUE)
FP_any <- rowSums((obs_state != -9) & (true_state == -9), na.rm = TRUE)
TN_any <- rowSums((obs_state == -9) & (true_state == -9), na.rm = TRUE)
FN_any <- rowSums((obs_state == -9) & (true_state != -9), na.rm = TRUE)


# CDR test performance for MCI vs HEALTHY
TP_mci <- rowSums((obs_state == 0) & (true_state == 0), na.rm = TRUE)
FP_mci <- rowSums((obs_state == 0) & (true_state == -9), na.rm = TRUE)
TN_mci <- rowSums((obs_state == -9) & (true_state == -9), na.rm = TRUE)
FN_mci <- rowSums((obs_state == -9) & (true_state == 0), na.rm = TRUE)

# CDR test performance for DEMENTIA vs MCI
TP_dem <- rowSums((obs_state > 0) & (true_state > 0), na.rm = TRUE)
FP_dem <- rowSums((obs_state > 0) & (true_state == 0), na.rm = TRUE)
TN_dem <- rowSums((obs_state == 0) & (true_state == 0), na.rm = TRUE)
FN_dem <- rowSums((obs_state == 0) & (true_state >= 0), na.rm = TRUE)


healthy_v_any <- test_performance(
  TP = TP_any,
  FP = FP_any,
  TN = TN_any,
  FN = FN_any
)

mci_v_healthy <- test_performance(
  TP = TP_mci,
  FP = FP_mci,
  TN = TN_mci,
  FN = FN_mci
)

dementia_v_mci <- test_performance(
  TP = TP_dem,
  FP = FP_dem,
  TN = TN_dem,
  FN = FN_dem
)
