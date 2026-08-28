## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("S1PCPPOS_EMR"),
  description = "Assess BHA eligibility via eRADAR each year until a PCP diagnosis (imperfect PCP)",
  
  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  # eRADAR-style EHR risk flag: assumed selection probabilities, set in this config.
  probs_select = f.select_matrix(h = 0.09, mci = 0.8, dem = 1),
  rr.select_prior = 2,    # RR of reporting concern again after reporting it last cycle
  repeat_interval = 1,    # years between BHA administrations
  prob_pcpfu = 1,
  
  # IMPERFECT PCP
  p.PCP_confirm_TP = c(0.75, 0.88, 0.95, 0.98), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 0.65,   # P(PCP correctly dismisses a BHA false positive)
  
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }  
  
)