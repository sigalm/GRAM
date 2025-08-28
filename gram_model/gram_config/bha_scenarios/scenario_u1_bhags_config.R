## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("U1_BHAGS"),
  description = "Repeat screening each year with the full BHA including informant survey",
  
  # Test parameters
  test = "BHA-GS",
  sensitivity = l.inputs[["sens_BHAGS"]], 
  specificity = l.inputs[["spec_BHAGS"]],
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 60,                   # age at which first BHA is administered
  probs_cogcon = l.inputs[["m.cogcon"]],  # Defaults to everyone getting tested. Use m.cogcon_spon or m.cogcon_elic for alternatives.
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    FALSE         # default: no stopping until death
  }  
  
)