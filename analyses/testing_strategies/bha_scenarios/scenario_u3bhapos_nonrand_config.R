## Scenario Config File

probs_select <-  l.inputs[["m.cogcon"]]
probs_select[,"h"] <- 0.315
probs_select[,"mci"] <- 0.63
probs_select[,"dem"] <- 0.63

scenario_inputs <- list(
  
  title       = as.factor("U3BHAPOS-NONRAND"),
  description = "Repeat BHA every 3 years for 1/3 of cohort with self-selection of CI, until first positive BHA",

  # Test parameters
  test        = "BHA-GS", 
  sensitivity = l.inputs[["sens_BHAGS"]], 
  specificity = l.inputs[["spec_BHAGS"]],
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,                     # age at which last BHA is administered
  probs_cogcon = probs_select,  # Defaults to everyone getting tested. Use m.cogcon_reactive or m.cogcon_selective for alternatives.
  repeat_interval = 3,    # years between BHA administrations
  prob_pcpfu = 0,         # Explicitly set to 0 to prevent default PCP follow-up behavior
  cohort_split = 3,       # 1/3 of the cohort gets tested every year (thus everyone gets tested once every 3 years) 
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         
  }  
  
)