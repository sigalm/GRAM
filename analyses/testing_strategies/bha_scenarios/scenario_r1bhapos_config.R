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
  # Spontaneous, patient-initiated concern. Carried over verbatim from the retired
  # data/cogcon/m.cogcon_reactive.RDS (KP preliminary data).
  probs_select = f.select_matrix(h = 0.01, mci = 0.1, dem = 0.3),
  rr.select_prior = 2,    # RR of reporting concern again after reporting it last cycle
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