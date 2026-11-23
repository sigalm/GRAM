## Scenario Config File

# Acceptance is the eRADAR arm's, built from the eRADAR flag and the prevalences it was
# derived in, so the two selective arms differ only in who is selected and in whether
# selection persists. See scenario_s1bhapos_emr_config.R.
eradar <- readRDS("data/eradar/eradar_flag_by_state.RDS")
eradar_flag <- f.select_matrix(h = eradar$h, mci = eradar$mci, dem = eradar$dem)

# Selection: P(endorses the question | state), 0.110 / 0.413 / 0.986, derived once by
# data/question/generate_question_endorse_by_state.R from the GRAMish workbook's P(tested)
# divided by the eRADAR take-up below; sources and method are there and in data/README.md.
question <- readRDS("data/question/question_endorse_by_state.RDS")

scenario_inputs <- list(

  title       = as.factor("S1BHAPOS-QUESTION"),
  description = "Assess BHA eligibility via a prompt each year until first positive result",

  # Test parameters
  test        = "BHA-GS",
  sensitivity = BHA_GS$sens,
  specificity = BHA_GS$spec,

  # Global parameters
  HCARE = 1,             # 1 = requires healthcare provider, 0 = ignore

  # Core scenario parameters
  age_first_test  = 65,                   # age at which first BHA is administered
  age_stop_test = 80,
  # Selective: question-based prompt at an annual wellness visit. The workbook's
  # 0.041 / 0.342 / 0.986 is P(endorses AND takes the test up); this is the endorsing part
  # alone, so that the take-up below is not applied twice. With it, P(tested) at a first
  # offer is the workbook's 0.041 / 0.342 / 0.986 again.
  # Workbook TCI folds into h; workbook memory loss folds into mci.
  probs_select = f.select_matrix(h = question$h, mci = question$mci, dem = question$dem),
  # Unlike the eRADAR flag, endorsing the question is something the patient does anew at
  # each visit, so selection does not persist: everyone due is redrawn at probs_select.
  select_persists = FALSE,
  # P(takes the test up | selected), by true state: the eRADAR arm's matrix, unchanged
  # (base case 0.373 / 0.828 / 1). It is not re-solved against this arm's selection, so
  # the share impaired among those tested is not held at 45% here: it is about 50% at a
  # first offer. Changing the uptake or p_impaired_tested here means regenerating the
  # endorsement file, which was divided by exactly this matrix.
  probs_accept = f.accept_matrix(eradar_flag,
                                 uptake            = 0.50,   # base case
                                 p_impaired_tested = 0.45,
                                 prev_mci          = eradar$prev_mci,
                                 prev_dem          = eradar$prev_dem),
  # Someone who declined their last offer takes the next one up at this, whatever their
  # state, as in the eRADAR arm. They are asked again the next year, but must endorse the
  # question again to be offered; the 0.2 applies whenever that next offer comes.
  # Source: email from Jim; citation to be added.
  p.accept_after_decline = 0.2,
  # A negative test means no re-test for 3 years, as in the eRADAR arm. A decline does
  # not start the interval, and neither does not being selected: anyone due and not
  # selected is asked again the next year.
  repeat_interval = 3,    # years between BHA administrations
  # No rr.select_prior: a prior test result does not change the next draw.

  # Optional parameters
  NP = NULL,

  # Stop parameters
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE         # default: no stopping until death
  }

)
