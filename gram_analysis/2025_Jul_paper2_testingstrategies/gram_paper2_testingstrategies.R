source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
source("gram_calibration/gram_benchmarking_helpers.r")
library(tableone)
sample1 <- readRDS("gram_data/acs_data/acs_age50.rds")



l.inputs1 <- l.inputs
l.inputs1[["n.cycle"]] <- 51
l.inputs1[["n.ind"]] <- 1000
l.inputs1[["p.HCARE_start"]] <- c(0.25,0.75)
l.inputs1[["hr.mort_mci_age"]] <- c(1, 1, 1)
l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(1, 1, 1)
l.inputs1[["seed_stochastic"]] <- 20250624

l.inputs1[["r.CDRslow_mean"]] <-  
  (seq(0, 1, length.out = 51)^1.5) * (1.5 * l.inputs[["r.CDRslow_mean"]])
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 2.2 /4


l.inputs1$scenario

universal <- f.wrap_run(l.inputs1, microdata = sample1)
