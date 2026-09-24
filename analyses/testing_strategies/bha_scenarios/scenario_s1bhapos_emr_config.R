## Scenario Config File

scenario_inputs <- list(
  
  title       = as.factor("S1BHAPOS_EMR"),
  description = "Assess BHA eligibility via eRADAR each year until first positive result",
  
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
  # Source: GRAMish workbook, "P conditional testing" tab, 14 Sep 2026.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = 0.041, mci = 0.371, dem = 0.734),
  # The flag belongs to the record: once selected, never redrawn. Whether an eRADAR
  # flag ever clears is conjecture the data do not settle; a sticky flag is the
  # working assumption for now.
  select_persists = TRUE,
  # A negative test does not clear the flag, but the person is not re-tested for 3
  # years, mirroring the inclusive arm. repeat_interval counts from the last TEST, so
  # a flagged person who did not take up the offer is offered again the next year, at
  # the unadjusted probability.
  repeat_interval = 3,    # years between BHA administrations
  # No rr.select_prior: repeats are independent. Neither a declined offer nor a prior
  # test result changes the next draw. probs_select is compound P(flagged) x
  # P(accept | flagged), so SELECT == 0 cannot separate "not flagged" from "flagged,
  # declined", and a decline penalty would hit both; left out until that split is
  # identified from the trial data.
  
  # Optional parameters
  NP = NULL,              
  
  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         # default: no stopping until death
  }  
  
)