## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("R1BHAPOS"),
  description = "Repeat BHA whenever concerns are brought up, until first positive result",

  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,                     # age at which to stop testing
  # Reactive: spontaneous, patient-initiated concern.
  # Source: GRAMish workbook, "P conditional testing" tab, 14 Sep 2026.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = 0.004, mci = 0.174, dem = 0.783),
  # No rr.select_prior: selection is redrawn at the unadjusted conditional probability at
  # every assessment. A prior test result carries no reassurance or persistence effect.
  repeat_interval = 1,    # years between BHA administrations
  
  # Optional parameters
  NP = NULL,              
  pause_after_FP = NULL,  
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE
  }  
  
)