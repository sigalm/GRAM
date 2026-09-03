## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("R1PCPPOS"),
  description = "Repeat BHA whenever concerns are brought up until a PCP diagnosis",
  
  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  # Reactive: spontaneous, patient-initiated concern.
  # Source: GRAMish workbook v25, 6 Feb 2026. Estimated prior to the July 2026
  # recalibration, so the undiagnosed case mix they were solved against has shifted.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = 0.005, mci = 0.2, dem = 0.9),
  # No rr.select_prior: selection is redrawn at the unadjusted conditional probability at
  # every assessment. A prior test result carries no reassurance or persistence effect.
  repeat_interval = 1,    # years between BHA administrations
  prob_pcpfu = 1,         # prob PCP follow up given a POSITIVE BHA
  # Previously inherited from model/setup.R; written out here now that strategy
  # parameters live only in the config. Values unchanged.
  p.PCP_confirm_TP = c(0.75, 0.88, 0.95, 0.98), # mci, mild, moderate, severe dem
  p.PCP_reject_FP  = 0.65,   # P(PCP correctly dismisses a BHA false positive)
  
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }  
  
)