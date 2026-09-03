## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("U3BHAPOS-RAND50"),
  description = "Repeat BHA every 3 years for 1/3 of cohort regardless of cognitive concern until first positive result",

  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,                     # age at which last BHA is administered
  # Inclusive: 50% random opt-in, independent of cognitive status.
  # Source: GRAMish workbook v25, 6 Feb 2026. Estimated prior to the July 2026
  # recalibration, so the undiagnosed case mix they were solved against has shifted.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = 0.5, mci = 0.5, dem = 0.5),
  # No rr.select_prior: selection is redrawn at the unadjusted conditional probability at
  # every assessment. A prior test result carries no reassurance or persistence effect.
  repeat_interval = 3,    # years between BHA administrations
  prob_pcpfu = 0,         # Explicitly set to 0 to prevent default PCP follow-up behavior
  cohort_split = 3,       # 1/3 of the cohort gets tested every year (thus everyone gets tested once every 3 years) 
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         
  }  
  
)