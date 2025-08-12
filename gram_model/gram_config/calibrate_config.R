# Calibration configuration

calibrate <- function(inputs, n) {
  
  inputs[["n.cycle"]] <- 51
  inputs[["n.ind"]] <- n
  inputs[["p.HCARE_start"]] <- c(0.25,0.75)
  inputs[["hr.mort_mci_age"]] <- c(1, 1, 1)
  inputs[["hr.mort_mod_age"]] <- inputs[["hr.mort_sev_age"]] <- c(1, 1, 1)
  inputs[["seed_stochastic"]] <- 20250624
  
  inputs[["r.CDRslow_mean"]] <-  
    (seq(0, 1, length.out = 51)^1.5) * (1.5 * inputs[["r.CDRslow_mean"]])
  inputs[["m.hr_mci"]] <- inputs[["m.hr_mci"]] * 2.2 /4
  
  return(inputs)
  
}
  


