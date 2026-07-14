# Calibration configuration

calibrate <- function(inputs, n) {
  
  inputs[["n.cycle"]] <- 51
  inputs[["n.ind"]] <- n
  inputs[["p.HCARE_start"]] <- c(0.25,0.75)
  inputs[["hr.mort_mci_age"]] <- c(1, 1, 1)
  inputs[["hr.mort_mod_age"]] <- inputs[["hr.mort_sev_age"]] <- c(1, 1, 1)
  inputs[["seed_stochastic"]] <- 20250624
  
  inputs[["param1"]]  <- 1.7900
  inputs[["param2a"]] <- 1.8000
  inputs[["param2b"]] <- 1.6680
  inputs[["m.hr_mci"]] <- inputs[["m.hr_mci"]] * inputs[["param1"]]
  inputs[["r.CDRslow_mean"]] <- (seq(0, 1, length.out = 51)^inputs[["param2a"]]) *
                                 (inputs[["param2b"]] * inputs[["r.CDRslow_mean"]])
  
  return(inputs)
  
}
  


