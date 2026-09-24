## Scenario Config File

# Selective: eRADAR-style EHR risk flag, P(flagged | state).
# Source: GRAMish workbook, "P conditional testing" tab, selective participation set to
# 100% so that P(tested) is P(flagged), 24 Sep 2026.
# Workbook TCI folds into h; workbook memory loss folds into mci.
eradar_flag <- f.select_matrix(h = 0.083, mci = 0.795, dem = 0.994)

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
  probs_select = eradar_flag,
  # The flag belongs to the record: once selected, never redrawn. Whether an eRADAR
  # flag ever clears is conjecture the data do not settle; a sticky flag is the
  # working assumption for now.
  select_persists = TRUE,
  # P(takes the test up | flagged), by true state, built from one overall uptake. See
  # f.accept_matrix() for how. The anchors and the MCI:dementia mix are from the
  # GRAMish workbook: a dementia:MCI testing ratio of 3.3 at the empirical 30% uptake,
  # chosen there so that P(tested | dem) matches the published eRADAR dementia
  # sensitivity, and background prevalence of 8.25% MCI and 0.93% dementia.
  probs_accept = f.accept_matrix(eradar_flag,
                                 uptake     = 0.50,   # base case
                                 ratio_ref  = 3.3,
                                 uptake_ref = 0.30,
                                 prev_mci   = 0.0825,
                                 prev_dem   = 0.0093),
  # A flagged person who declined their last offer takes the next one up at this,
  # whatever their state. Source: email from Jim; citation to be added.
  p.accept_after_decline = 0.2,
  # A negative test does not clear the flag, but the person is not re-tested for 3
  # years, mirroring the inclusive arm. A decline restarts the interval too, so a
  # flagged person who did not take up the offer is offered again 3 years later.
  repeat_interval = 3,    # years between BHA administrations
  # No rr.select_prior: a prior test result does not change the next draw, and with
  # select_persists the flagged are never redrawn anyway.
  prob_pcpfu = 1,
  
  # IMPERFECT PCP
  p.PCP_confirm_TP = c(0.75, 0.88, 0.95, 0.98), # mci, mild, moderate, severe dem
  p.PCP_reject_FP = 0.65,   # P(PCP correctly dismisses a BHA false positive)
  
  
  # Optional parameters
  NP = NULL,              
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_PCP_pos == 1           # do not repeat BHA if positive PCP assessment
  }  
  
)