## Scenario Config File - TEMPLATE
## Copy this file and rename to: scenario_<NAME>_config.R

scenario_inputs <- list(
  
  # ============================================================================
  # SCENARIO IDENTIFICATION
  # ============================================================================
  title       = as.factor("SCENARIO_NAME"),
  description = "Brief description of the testing scenario",

  # ============================================================================
  # TEST PARAMETERS
  # ============================================================================
  test        = "BHA-GS",                       # Test type
  sensitivity = BHA_GS$sens,       # Test sensitivity
  specificity = BHA_GS$spec,       # Test specificity
  
  # ============================================================================
  # GLOBAL PARAMETERS
  # ============================================================================
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore
  
  # ============================================================================
  # CORE SCENARIO PARAMETERS
  # ============================================================================
  age_first_test  = 65,                         # Age at which first BHA is administered
  age_stop_test   = 80,                         # Age at which to stop testing
  # Universal: everyone assessable is tested; selection is not a gate.
  probs_select = f.select_matrix(h = 1, mci = 1, dem = 1),
                                                    # or custom matrix with dimensions age x cog status
  # Selection is redrawn at the unadjusted conditional probability at every assessment
  # unless one of the next two is set. Both are optional and both default to off.
  # select_persists = TRUE,        # once selected, never redrawn -- for a flag that
                                   # belongs to the record (eRADAR) rather than to a
                                   # decision the patient remakes each year
  # rr.select_prior = c(neg = 0.5, pos = 1),
                                   # RR on re-selection, keyed on the result of a test in
                                   # the IMMEDIATELY PRECEDING cycle. Names are optional
                                   # (neg, pos); an omitted name means no effect for that
                                   # result. No test last cycle is always an unadjusted
                                   # redraw, so this is inert wherever repeat_interval > 1.
  prob_pcpfu      = NULL,                       # Probability of PCP follow-up after cognitive test (default is no follow-up)
  repeat_interval = 1,                          # Years between BHA administrations
  cohort_split    = NULL,                       # Split cohort into groups (useful for alternated testing scenarios, e.g., test half the cohort every other year)
  
  # ============================================================================
  # OPTIONAL PARAMETERS
  # ============================================================================
  NP              = NULL,                       # Neuropsychological testing
  pause_after_FP  = NULL,                       # Pause testing after false positive (default is no test after positive)
  
  # ============================================================================
  # STOP RULE
  # ============================================================================
  # Available args in stop_rule function:
  #   any_BHA_pos    - Any prior positive BHA result
  #   NP             - Neuropsychological test result (lagged)
  #   PCP.lag        - PCP evaluation result (lagged)
  #   repeat_after_FP - Repeat after false positive flag
  
  stop_rule = function(...) {
    args <- list(...)
    
    # Example: Stop testing after first positive BHA
    args$any_BHA_pos == TRUE
  }
  
)
