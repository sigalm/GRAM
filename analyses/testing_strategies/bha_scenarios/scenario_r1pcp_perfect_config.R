## Scenario Config File

scenario_inputs <- list(

  title       = as.factor("R1PCPPOS_PERFECT"),
  description = "Repeat BHA whenever concerns are brought up until a PCP diagnosis (perfect PCP)",

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
  # Everyone who raises a concern is tested: the concern is the patient's own
  # initiative, so there is no offer to decline. Set explicitly so every arm states its
  # uptake.
  probs_accept = f.select_matrix(h = 1, mci = 1, dem = 1),
  # No rr.select_prior: selection is redrawn at the unadjusted conditional probability at
  # every assessment. A prior test result carries no reassurance or persistence effect.
  repeat_interval = 1,    # years between BHA administrations
  prob_pcpfu = 1,

  # PERFECT PCP
  p.PCP_confirm_TP = c(1, 1, 1, 1), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 1,   # P(PCP correctly dismisses a BHA false positive)


  # Optional parameters
  NP = NULL,

  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }

)
