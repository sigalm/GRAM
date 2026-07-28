######################################## GRAM HELPER FUNCTIONS: EPIDEMIOLOGY FUNCTIONS ########################################
# This section defines the functions for epidemiologic variables.


######################################## MCI PROBABILITY

f.calc_MCIprob <- function(l.inputs, v.AGE, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag) {
  # As in the life table, hold the incidence lookup at the oldest tabulated age so that cohorts
  # which age beyond it (or start above 50) stay in bounds.
  age_index <- f.age_index(v.AGE, nrow(l.inputs[["m.hr_mci"]]))
  hazards.age <- l.inputs[["m.hr_mci"]][matrix(data = age_index, ncol = 1)]
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
