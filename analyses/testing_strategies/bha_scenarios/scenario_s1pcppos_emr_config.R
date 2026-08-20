## Scenario Config File

probs_select <-  l.inputs[["m.cogcon"]]
probs_select[,"h"] <- 0.09
probs_select[,"mci"] <- 0.8
probs_select[,"dem"] <- 1

scenario_inputs <- list(
  
  title       = as.factor("S1PCPPOS_EMR"),
  description = "Assess BHA eligibility via eRADAR each year until a PCP diagnosis (imperfect PCP)",
  
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
  prob_pcpfu = 1,
  
  # IMPERFECT PCP
  p.PCP_confirm_TP = c(0.75, 0.88, 0.95, 0.98), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 0.65,   # P(PCP correctly dismisses a BHA false positive)
  
  
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