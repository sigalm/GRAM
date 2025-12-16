## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("S1PCPPOS"),
  description = "Assess BHA eligibility with a prompt each year until a PCP diagnosis",
  
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
  probs_cogcon = l.inputs[["m.cogcon_elic"]],  # Defaults to everyone getting tested. Use m.cogcon_spon or m.cogcon_elic for alternatives.
  repeat_interval = 1,    # years between BHA administrations
  prob_pcpfu = 1,
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }  
  
)