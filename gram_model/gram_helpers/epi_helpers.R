######################################## GRAM HELPER FUNCTIONS: EPIDEMIOLOGY FUNCTIONS ########################################
# This section defines the functions for epidemiologic variables.


######################################## MCI PROBABILITY

f.calc_MCIprob <- function(l.inputs, v.AGE, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag) {
  hazards.age <- l.inputs[["m.hr_mci"]][matrix(data = round(v.AGE,0), ncol = 1) - 50 + 1]
  hazard <- hazards.age * exp(
    l.inputs[["log_EDU"]] * v.EDU.lag +
      l.inputs[["log_SEX"]] * v.SEX.lag +
      l.inputs[["log_RACEETHblack"]] * (v.RACEETH.lag == 1) +
      l.inputs[["log_RACEETHhisp"]] * (v.RACEETH.lag == 2) +
      l.inputs[["log_APOE4"]] * v.APOE4.lag +
      l.inputs[["log_MEDBUR"]] * v.MEDBUR.lag +
      l.inputs[["log_INCOMEmed"]] * (v.INCOME.lag == 1) +
      l.inputs[["log_INCOMEhi"]] * (v.INCOME.lag == 2))
  
  prob_mci <- (1 - exp(-hazard)) * l.inputs[["rr.Px_mci"]]
  
  return(prob_mci)
}
