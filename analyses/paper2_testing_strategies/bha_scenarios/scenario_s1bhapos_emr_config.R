## Scenario Config File

probs_select <-  l.inputs[["m.cogcon"]]
probs_select[,"h"] <- 0.09
probs_select[,"mci"] <- 0.8
probs_select[,"dem"] <- 1



scenario_inputs <- list(
  
  title       = as.factor("S1BHAPOS_EMR"),
  description = "Assess BHA eligibility via eRADAR each year until first positive result",
  
  # Test parameters
  test        = "BHA-GS", 
  sensitivity = l.inputs[["sens_BHAGS"]], 
  specificity = l.inputs[["spec_BHAGS"]],
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  probs_cogcon = probs_select,
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         # default: no stopping until death
  }  
  
)