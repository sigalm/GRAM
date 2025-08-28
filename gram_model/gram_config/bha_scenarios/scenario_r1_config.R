## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("R1"),
  description = "Repeat BHA when concerns are brought up regardless of test history",
  
  # Test parameters
  test        = "BHA-CS", 
  sensitivity = l.inputs[["sens_BHACS"]], 
  specificity = l.inputs[["spec_BHACS"]],
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 60,                   # age at which first BHA is administered
  probs_cogcon = l.inputs[["m.cogcon_spon"]],  # Defaults to everyone getting tested. Use m.cogcon_spon or m.cogcon_elic for alternatives.
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    FALSE == TRUE         # default: no stopping until death
  }  
  
)