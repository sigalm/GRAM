## Scenario Config File


scenario_inputs <- list(
  
  title       = as.factor("U1-HYBRID_FU"),
  description = "Assess BHA eligibility via a prompt, with some % of patients receive PCP follow-up",
  
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
  prob_pcpfu = 0.2,  # the probability that a patient with a regular healthcare provider will receive PCP follow up after a cognitive test
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  PET = NULL,             
  pause_after_FP = NULL,  
  
  # Stop parameters
  
  
  # default stop_rule args:
  # any_BHA_pos = v.any_BHA_pos, 
  # NP = v.NP.lag, 
  # PET = v.PET.lag,
  # PCP.lag = v.PCP.lag,   
  # repeat_after_FP = scenario$repeat_after_FP
  
  stop_rule = function(...) {
    args <- list(...)
    
    stop_test <- args$any_BHA_pos == TRUE    # default action is to look at BHA test result
    stop_test[!is.na(args$PCP.lag)] <- (args$PCP.lag[!is.na(args$PCP.lag)] == 1) # override BHA test result if PCP evaluation was done
    
    return(stop_test)
  }
  
)