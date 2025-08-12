## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("U1"),
  description = "Annual universal screening regardless of cognitive concerns or test history",
  test        = TRUE,                 # BHA test is on

  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 60,                   # age at which first BHA is administered
  probs_cogcon = l.inputs[["m.cogcon"]],  # Defaults to everyone getting tested. Use m.cogcon_spon or m.cogcon_elic for alternative scenarios.
  
  # Initiation parameters
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Test interval parameters
  repeat_interval = 1,    # years between BHA administrations
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    FALSE == TRUE         # evaluates to FALSE every time, so no stop until death
  }  
  
)