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
  # Inclusive: everyone due is offered a test, so selection is not a gate; who is tested
  # is decided by uptake, 50% random opt-in, independent of cognitive status.
  probs_select = f.select_matrix(h = 1, mci = 1, dem = 1),
  # P(takes the test up | offered), by true state.
  # Source: GRAMish workbook v25, 6 Feb 2026. Estimated prior to the July 2026
  # recalibration, so the undiagnosed case mix they were solved against has shifted.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_accept = f.select_matrix(h = 0.5, mci = 0.5, dem = 0.5),
  # Someone who declined their last offer takes the next one up at this, whatever their
  # state. Source: email from Jim; citation to be added. A decline restarts
  # repeat_interval, so the next offer is 3 years on -- when cohort_split would have
  # made them due again anyway.
  p.accept_after_decline = 0.2,
  # No rr.select_prior: a prior test result carries no reassurance or persistence effect.
  repeat_interval = 3,    # years between BHA administrations
  prob_pcpfu = 0,         # Explicitly set to 0 to prevent default PCP follow-up behavior
  cohort_split = 3,       # 1/3 of the cohort gets tested every year (thus everyone gets tested once every 3 years) 
  
  # Optional parameters
  NP = NULL,              
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         
  }  
  
)