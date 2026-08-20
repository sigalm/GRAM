######################################## GRAM HELPER FUNCTIONS: EPIDEMIOLOGY FUNCTIONS ########################################
# This section defines the functions for epidemiologic variables.


######################################## MCI PROBABILITY

f.calc_MCIprob <- function(l.inputs, v.AGE, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag) {
  # As in the life table, hold the incidence lookup at the oldest tabulated age so that cohorts
  # which age beyond it (or start above 50) stay in bounds.
  # param1 is the calibrated multiplier on the published age-specific incidence hazards. It is
  # applied here, at the point of use, rather than being written back over m.hr_mci in l.inputs:
  # that older arrangement destroyed the published values and was not idempotent, so applying it
  # twice silently squared the multiplier.
  age_index <- f.age_index(v.AGE, nrow(l.inputs[["m.hr_mci"]]))
  hazards.age <- l.inputs[["m.hr_mci"]][matrix(data = age_index, ncol = 1)] * l.inputs[["param1"]]
  hazard <- hazards.age * exp(
    l.inputs[["log_EDU"]] * v.EDU.lag +
      l.inputs[["log_SEX"]] * (v.SEX.lag == 2) +
      l.inputs[["log_RACEETHblack"]] * (v.RACEETH.lag == 1) +
      l.inputs[["log_RACEETHhisp"]] * (v.RACEETH.lag == 2) +
      l.inputs[["log_APOE4"]] * v.APOE4.lag +
      l.inputs[["log_MEDBUR"]] * v.MEDBUR.lag +
      l.inputs[["log_INCOMEmed"]] * (v.INCOME.lag == 1) +
      l.inputs[["log_INCOMEhi"]] * (v.INCOME.lag == 2))
  
  prob_mci <- (1 - exp(-hazard)) * l.inputs[["rr.Px_mci"]]

  return(prob_mci)
}


######################################## PERSON-LEVEL TIME TO EVENT

# Collapse a.out into one row per individual so that it can be used for survival analysis (e.g., epiR). 
# Include time to first event, and what ended their follow-up.

#
#   event = "dementia"  first cycle with SEV >= 1
#   event = "mci"       first cycle with SYN == 1, i.e. onset of impairment. Anyone whose
#                       first impaired cycle is already at dementia level still counts, so
#                       no onset is missed if an individual skips past the MCI band.
#
# Returns baseline attributes (all of v.attr_names, as at baseline_cycle) plus:
#   id      column index in a.out
#   time    years from baseline to whatever ended follow-up
#   status  0 = censored alive and event-free at the end of the window
#           1 = event
#           2 = died event-free (a competing risk, not a censoring)
#
# Individuals already in the event state at baseline, or not alive at baseline, are dropped:
# they were never at risk.
#
# Times are years. The model only observes state at cycle boundaries, so an event in
# cycle t is recorded at t - baseline_cycle, the end of the interval it happened in.

f.person_time <- function(a.out, event = c("dementia", "mci"), baseline_cycle = 1,
                          max_followup = NULL) {

  event <- match.arg(event)

  last_cycle <- dim(a.out)[1]
  if (!is.null(max_followup)) last_cycle <- min(last_cycle, baseline_cycle + max_followup)
  if (last_cycle <= baseline_cycle) stop("no follow-up cycles between baseline_cycle and max_followup")

  cycles <- (baseline_cycle + 1):last_cycle
  n.fu   <- length(cycles)

  # matrix() rather than plain [ , , ] so a single follow-up cycle still gives a matrix
  m.SYN   <- matrix(a.out[cycles, "SYN", ],   nrow = n.fu)
  m.SEV   <- matrix(a.out[cycles, "SEV", ],   nrow = n.fu)
  m.ALIVE <- matrix(a.out[cycles, "ALIVE", ], nrow = n.fu)

  # SEV is NA unless SYN == 1, and all attributes but ALIVE are NA once someone has died,
  # so NA here means "not a case this cycle" either way.
  m.event <- if (event == "dementia") m.SEV >= 1 else m.SYN == 1
  m.event[is.na(m.event)] <- FALSE

  f.first_row <- function(m) apply(m, 2, function(x) if (any(x)) which(x)[1] else NA_integer_)

  row_event <- f.first_row(m.event)
  row_death <- f.first_row(m.ALIVE == 0)

  # Mortality is updated before health within a cycle, so someone who dies in cycle t has SEV
  # NA there and cannot also be an incident case.
  status <- ifelse(!is.na(row_event) & (is.na(row_death) | row_event <= row_death), 1L,
                   ifelse(!is.na(row_death), 2L, 0L))
  time   <- ifelse(status == 1L, row_event,
                   ifelse(status == 2L, row_death, n.fu))

  at_baseline <- if (event == "dementia") {
    a.out[baseline_cycle, "SEV", ] >= 1
  } else {
    a.out[baseline_cycle, "SYN", ] == 1
  }
  at_baseline[is.na(at_baseline)] <- FALSE
  at_risk <- which(a.out[baseline_cycle, "ALIVE", ] == 1 & !at_baseline)

  out <- as.data.frame(t(a.out[baseline_cycle, , at_risk]))
  out$id     <- at_risk
  out$time   <- time[at_risk]
  out$status <- status[at_risk]

  return(out)
}
