## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("S1PCPPOS_PERFECT_EMR"),
  description = "Assess BHA eligibility via eRADAR each year until a PCP diagnosis (perfect PCP)",
  
  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  # Selective: eRADAR-style EHR risk flag.
  # Source: GRAMish workbook v25, 6 Feb 2026. Estimated prior to the July 2026
  # recalibration, so the undiagnosed case mix they were solved against has shifted.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = 0.09, mci = 0.8, dem = 1),
  # Eligibility is reassessed from scratch every cycle: no select_persists, no
  # rr.select_prior. Whether an eRADAR flag should stick, and whether a negative result
  # should reassure, are both conjecture the data do not settle. A clean redraw is the
  # assumption that needs no defending, and the pathways it leaves out largely cancel:
  # those who test positive drop out of the pool anyway, and for those who test negative
  # a cleared flag IS a pure redraw.
  repeat_interval = 1,    # years between BHA administrations
  prob_pcpfu = 1,
  
  # PERFECT PCP
  p.PCP_confirm_TP = c(1, 1, 1, 1), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 1,   # P(PCP correctly dismisses a BHA false positive)
  
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }  
  
)