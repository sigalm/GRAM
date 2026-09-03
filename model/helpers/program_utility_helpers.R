######################################## GRAM HELPER FUNCTIONS: PROGRAM UTILITY HELPERS ########################################
# This script defines program utility helpers.

# Build a probs_select matrix.
# This is the ONLY place the 50:100 age grid is written down. f.update_SELECT() looks the
# probability up by row POSITION (round(AGE) - 50 + 1), not by the value in the age column,
# so the grid has to stay pinned to the model's age origin.
#
# h covers everyone with SYN < 1, which includes TCI; mci is SEV == 0; dem is SEV >= 1.
# Pass length-51 vectors instead of scalars for an age-varying probability.
f.select_matrix <- function(h, mci, dem, ages = 50:100) {
  data.frame(age = ages, h = h, mci = mci, dem = dem)
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

  m <- scen[["probs_select"]]
  if (!all(c("age", "h", "mci", "dem") %in% names(m)) || nrow(m) != 51) {
    stop(label, ": probs_select must have columns age/h/mci/dem and 51 rows (ages 50-100). ",
         "Build it with f.select_matrix().", call. = FALSE)
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