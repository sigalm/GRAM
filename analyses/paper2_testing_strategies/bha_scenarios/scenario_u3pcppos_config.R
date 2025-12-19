## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("U1PCPPOS"),
  description = "Repeat BHA every year regardless of cognitive concern until a PCP diagnosis",

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
  probs_cogcon = l.inputs[["m.cogcon"]],  # Defaults to everyone getting tested. Use m.cogcon_reactive or m.cogcon_selective for alternatives.
  prob_pcpfu = 1,  # the probability that a patient with a regular healthcare provider will receive PCP follow up after a cognitive test
  repeat_interval = 3,    # years between BHA administrations
  cohort_split = 3,
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1          # do not repeat BHA if positive PCP assessment
  }  
  
)