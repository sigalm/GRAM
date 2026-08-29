## Scenario Config File


scenario_inputs <- list(
  
  title       = as.factor("U1-HYBRID_FU"),
  description = "Assess BHA eligibility via a prompt, with some % of patients receive PCP follow-up",
  
  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  # Universal: everyone assessable is tested; selection is not a gate.
  probs_select = f.select_matrix(h = 1, mci = 1, dem = 1),
  rr.select_prior = 2,    # RR of reporting concern again after reporting it last cycle
  prob_pcpfu = 0.2,  # the probability that a patient with a regular healthcare provider will receive PCP follow up after a cognitive test
  # Previously inherited from model/setup.R; written out here now that strategy
  # parameters live only in the config. Values unchanged.
  p.PCP_confirm_TP = c(0.75, 0.88, 0.95, 0.98), # mci, mild, moderate, severe dem
  p.PCP_reject_FP  = 0.65,   # P(PCP correctly dismisses a BHA false positive)
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  
  
  # default stop_rule args:
  # any_BHA_pos = v.any_BHA_pos, 
  # NP = v.NP.lag, 
  # PCP.lag = v.PCP.lag,   
  # repeat_after_FP = scenario$repeat_after_FP
  
  stop_rule = function(...) {
    args <- list(...)
    
    stop_test <- args$any_BHA_pos == TRUE    # default action is to look at BHA test result
    stop_test[!is.na(args$PCP.lag)] <- (args$PCP.lag[!is.na(args$PCP.lag)] == 1) # override BHA test result if PCP evaluation was done
    
    return(stop_test)
  }
  
)