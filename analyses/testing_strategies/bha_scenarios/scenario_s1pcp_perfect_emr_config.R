## Scenario Config File

# Selective: eRADAR EHR risk flag, P(flagged | state), 0.122 / 0.404 / 0.659, derived
# once from the published eRADAR validation by data/eradar/generate_eradar_flag_by_state.R;
# sources and method are there and in data/README.md. The file also carries the
# prevalences the derivation used, which the acceptance below is built in too.
eradar <- readRDS("data/eradar/eradar_flag_by_state.RDS")
eradar_flag <- f.select_matrix(h = eradar$h, mci = eradar$mci, dem = eradar$dem)

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
  probs_select = eradar_flag,
  # The flag belongs to the record: once selected, never redrawn. Whether an eRADAR
  # flag ever clears is conjecture the data do not settle; a sticky flag is the
  # working assumption for now. The never-flagged are redrawn independently each year;
  # a sensitivity analysis correlating redraws with history is planned.
  select_persists = TRUE,
  # P(takes the test up | flagged), by true state: the same three assumptions as
  # scenario_s1bhapos_emr_config.R, where the reasoning is written out.
  probs_accept = f.accept_matrix(eradar_flag,
                                 uptake            = 0.50,   # base case
                                 p_impaired_tested = 0.45,
                                 prev_mci          = eradar$prev_mci,
                                 prev_dem          = eradar$prev_dem),
  # A flagged person who declined their last offer takes the next one up at this,
  # whatever their state -- dementia too, in place of its first-offer acceptance of 1.
  # Since no one with dementia declines a first offer, this reaches dementia only through
  # someone who declined while healthy or MCI and has since progressed.
  # Source: uptake among earlier non-participants in later rounds of colorectal cancer
  # screening: Kooyker et al., Int J Cancer 2020;147:1098-1106; Lo et al., Gut
  # 2015;64:282-291.
  p.accept_after_decline = 0.2,
  # A negative test does not clear the flag, but the person is not re-tested for 3
  # years, mirroring the inclusive arm. A decline does not start the interval: a
  # flagged person who did not take up the offer is offered again the next year.
  repeat_interval = 3,    # years between BHA administrations
  # No rr.select_prior: a prior test result does not change the next draw, and with
  # select_persists the flagged are never redrawn anyway.
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