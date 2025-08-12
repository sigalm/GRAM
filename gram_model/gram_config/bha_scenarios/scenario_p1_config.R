## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("P1"),
  description = "Assess BHA eligibility with a prompt each year, regardless of test history",
  test        = TRUE,                 # BHA test is on
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  DX    = 0,             # 1 = ignore, 0 = require no prior diagnosis
  
  # Core scenario parameters
  age_first_test  = 60,                   # age at which first BHA is administered
  probs_cogcon = l.inputs[["m.cogcon_elic"]],  # Defaults to everyone getting tested. Use m.cogcon_spon or m.cogcon_elic for alternatives.
  
  # Initiation parameters
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    FALSE == TRUE       # default: evaluates to FALSE, no stopping until death
  }  
  
)