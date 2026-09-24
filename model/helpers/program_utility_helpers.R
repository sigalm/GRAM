######################################## GRAM HELPER FUNCTIONS: PROGRAM UTILITY HELPERS ########################################
# This script defines program utility helpers.

# Build a probs_select or probs_accept matrix.
# This is the ONLY place the 50:100 age grid is written down. f.lookup_by_state() looks the
# probability up by row POSITION (round(AGE) - 50 + 1), not by the value in the age column,
# so the grid has to stay pinned to the model's age origin.
#
# h covers everyone with SYN < 1, which includes TCI; mci is SEV == 0; dem is SEV >= 1.
# Pass length-51 vectors instead of scalars for an age-varying probability.
f.select_matrix <- function(h, mci, dem, ages = 50:100) {
  data.frame(age = ages, h = h, mci = mci, dem = dem)
}


# Build a probs_accept matrix from an overall uptake, so that trying another uptake is
# one number in the config rather than another pass through the workbook.
#
# probs_select is P(flagged | state): the uptake = 100% column of the workbook. uptake is
# the share of the flagged who take the test up. Acceptance is P(test | flagged, state).
#
#   healthy    accepts at uptake. With the prevalence-weighted total fixed at
#              uptake x P(flagged), the impaired then average uptake too.
#   MCI / dem  split that average so P(test | dem) / P(test | MCI) = R(uptake), a power
#              curve through two anchors: ratio_ref at uptake_ref (the empirical
#              anchor) and the flag ratio dem / MCI at uptake = 100%, where every
#              acceptance is 1. R falls as uptake rises -- the dementia cases are mostly
#              in already, so more uptake brings in relatively more MCI -- and reaches
#              the flag ratio without a floor or a kink.
#
# prev_mci and prev_dem only enter as their ratio: the mix of MCI and dementia among
# the impaired. Every argument is required; the values belong in the scenario config.
f.accept_matrix <- function(probs_select, uptake, ratio_ref, uptake_ref, prev_mci, prev_dem) {
  stopifnot(uptake > 0, uptake <= 1, uptake_ref > 0, uptake_ref < 1)

  f_h <- probs_select$h; f_mci <- probs_select$mci; f_dem <- probs_select$dem

  ratio_100 <- f_dem / f_mci
  ratio <- ratio_100 * (ratio_ref / ratio_100)^(log(uptake) / log(uptake_ref))

  p_test_mci <- uptake * (f_mci * prev_mci + f_dem * prev_dem) / (prev_mci + ratio * prev_dem)
  p_test_dem <- ratio * p_test_mci

  m <- f.select_matrix(h   = rep(uptake, length(f_h)),
                       mci = p_test_mci / f_mci,
                       dem = p_test_dem / f_dem,
                       ages = probs_select$age)

  # Acceptance above 1 means these anchors ask the dementia flag for more tests than
  # there are flagged people with dementia at this uptake.
  if (any(m[, c("h", "mci", "dem")] > 1 + 1e-9)) {
    stop("f.accept_matrix: acceptance above 1 at uptake = ", uptake,
         ". Check ratio_ref / uptake_ref against probs_select.", call. = FALSE)
  }
  m
}


# Look up each person's probability in a matrix built by f.select_matrix(), by age and
# true state. Shared by selection and acceptance so the two can never disagree about
# which column a person falls in.
f.lookup_by_state <- function(m, v.AGE, v.SYN, v.SEV) {
  col <- case_when(
    v.SYN < 1   ~ 2,
    v.SEV == 0  ~ 3,
    v.SEV >= 1  ~ 4
  )
  m[matrix(data = c(round(v.AGE, 0) - 50 + 1, col), ncol = 2)]
}


# Map the result of a test in the previous cycle onto a relative risk of being selected
# again. rr.select_prior may be:
#   NULL                    no effect anywhere
#   a single unnamed number one effect after any prior result, positive or negative --
#                           what the pre-September-2026 configs meant by rr.select_prior = 2
#   a named vector          any of neg / pos, e.g. c(neg = 0.5)
# A name left out of a named vector means NO effect for that result, which is not the
# same as omitting the parameter: an arm can reassure after a negative and still leave a
# positive at the unadjusted probability. v.BHA.lag is -9 (not assessed last cycle) or
# -8 (assessed but not tested); neither is a result, so both come back at RR 1.
f.rr_by_prior_result <- function(rr.select_prior, v.BHA.lag) {

  rr <- rep(1, length(v.BHA.lag))
  if (is.null(rr.select_prior)) return(rr)

  if (is.null(names(rr.select_prior))) {
    rr[v.BHA.lag >= 0] <- rr.select_prior[[1]]
    return(rr)
  }

  result_code <- c(neg = 0, pos = 1)
  for (nm in names(rr.select_prior)) rr[v.BHA.lag == result_code[[nm]]] <- rr.select_prior[[nm]]
  rr
}


# Validate a scenario before it is run.
# A testing scenario has to carry its own parameters: there are deliberately no defaults
# in model/setup.R to fall back on, so a missing field is an error here rather than a
# confusing NULL subscript several minutes into f.run().
f.validate_scenario <- function(scen, label = "scenario") {
  if (is.null(scen[["test"]])) return(invisible(scen))   # natural history, nothing to check

  required <- c("probs_select", "sensitivity", "specificity")   # rr.select_prior is optional
  missing  <- required[vapply(required, function(f) is.null(scen[[f]]), logical(1))]
  if (length(missing)) {
    stop(label, " defines a test but is missing: ", paste(missing, collapse = ", "),
         ". Testing parameters must be set in the scenario config.", call. = FALSE)
  }

  for (nm in c("probs_select", "probs_accept")) {   # probs_accept is optional
    m <- scen[[nm]]
    if (is.null(m)) next
    if (!all(c("age", "h", "mci", "dem") %in% names(m)) || nrow(m) != 51) {
      stop(label, ": ", nm, " must have columns age/h/mci/dem and 51 rows (ages 50-100). ",
           "Build it with f.select_matrix().", call. = FALSE)
    }
  }

  # The after-decline probability only means something where offers can be declined.
  if (!is.null(scen[["p.accept_after_decline"]]) && is.null(scen[["probs_accept"]])) {
    stop(label, " sets p.accept_after_decline but not probs_accept. Without probs_accept ",
         "every offer is taken up, so there is never a decline for it to act on.", call. = FALSE)
  }

  # rr.select_prior is optional, but a name f.rr_by_prior_result does not recognise would
  # be silently ignored for the whole run, so it is caught here instead.
  rr <- scen[["rr.select_prior"]]
  if (!is.null(rr) && !is.null(names(rr))) {
    bad <- setdiff(names(rr), c("neg", "pos"))
    if (length(bad)) {
      stop(label, ": rr.select_prior has unrecognised name(s) ", paste(bad, collapse = ", "),
           ". Use neg and/or pos, keyed on the previous cycle's BHA result.", call. = FALSE)
    }
  }

  # PCP follow-up is optional, but asking for it without its probabilities is not.
  pcpfu <- scen[["prob_pcpfu"]]
  if (!is.null(pcpfu) && pcpfu > 0) {
    pcp_required <- c("p.PCP_confirm_TP", "p.PCP_reject_FP")
    pcp_missing  <- pcp_required[vapply(pcp_required, function(f) is.null(scen[[f]]), logical(1))]
    if (length(pcp_missing)) {
      stop(label, " sets prob_pcpfu = ", pcpfu, " but is missing: ",
           paste(pcp_missing, collapse = ", "), ".", call. = FALSE)
    }
  }

  invisible(scen)
}


# Load configuration files
load_scenario <- function(config_file, base_inputs) {
  env <- new.env()
  source(config_file, local = env)
  base_inputs[["scenario"]] <- modifyList(base_inputs[["scenario"]], env$scenario_inputs)
  f.validate_scenario(base_inputs[["scenario"]], label = basename(config_file))
  base_inputs
}


# Read the results of last saved run
latest_rds <- function(prefix, dir = output_dir) {
  files <- list.files(output_dir, pattern = paste0("scenario_", prefix, "_sim_.*\\.rds$"), full.names = TRUE)
  
  if (length(files) == 0) stop("No files found for prefix: ", prefix)
  latest <- sort(files, decreasing = TRUE)[1]
  
  message("Loading: ", latest)
  readRDS(latest)
}