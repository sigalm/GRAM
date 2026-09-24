## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("U1PCPPOS-PERFECT-RAND50"),
  description = "Repeat BHA every year regardless of cognitive concern until a PCP diagnosis (perfect PCP)",

  # Test parameters
  test        = "BHA-GS", 
  sensitivity = BHA_GS$sens, 
  specificity = BHA_GS$spec,
  
  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
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
  prob_pcpfu = 1,  # the probability that a patient with a regular healthcare provider will receive PCP follow up after a cognitive test
  
  # PERFECT PCP
  p.PCP_confirm_TP = c(1, 1, 1, 1), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 1,   # P(PCP correctly dismisses a BHA false positive)
  
  repeat_interval = 3,    # years between BHA administrations
  cohort_split = 3,
  
  # Optional parameters
  NP = NULL,              
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1          # do not repeat BHA if positive PCP assessment
  }  
  
)