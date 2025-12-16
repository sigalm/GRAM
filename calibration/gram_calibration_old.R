# Calibrate average rate of change in CDR-SB for slow progressors to achieve expected time in MCI
# 1. Try optim() with method "L-BFGS-B" -- focusing only on one parameter (slow mean) and not fully exploring the parameter space
# 2. Try nlminb() which is supposed to work better when there are multiple parameters -- same issue as above
# 3. Try GenSA (simulated annealing) to avoid local minima

library(GenSA)

l.inputs_optim <- l.inputs
l.inputs_optim[["r.CDRfast_mean"]] <- 1.00
l.inputs_optim[["r.CDRfast_sd1"]] <- 0
l.inputs_optim[["r.CDRslow_mean"]] <- 0.51
l.inputs_optim[["r.CDRslow_sd1"]] <- 0

target_mci_dur <- 2
target_dem_dur <- 4 
of <- function(params) {
  print(params)
  l.inputs_optim[["r.CDRslow_mean"]] <- params[1]
  # l.inputs_optim[["r.CDRfast_mean"]] <- params[2]
  result <- f.wrap_run(a.random,l.inputs, printLevel = 0)
  mci_dur <- result$aggregated_results$MCI.sum
  dem_dur <- sum(
    result$aggregated_results$SEV1.sum + result$aggregated_results$SEV2.sum + result$aggregated_results$SEV3.sum
  )
  mci_difference <- (mci_dur - target_mci_dur)^2
  dem_difference <- (dem_dur - target_dem_dur)^2
  tot_difference <- 0.7*mci_difference + 0.3*dem_difference
  return(mci_difference)
}

initial_guess <- c(1.0, 1.5)

result_optim <- optim(par=1.0, fn = of, method = "L-BFGS-B",
                      lower = 0.5, upper = 2)
best_optim <- result_optim$par
best_optim


result_gensa <- GenSA(par = initial_guess, fn = of,
                lower = c(0.5, 0.5), upper = c(2, 3),
                control = list(verbose = TRUE, maxit = 1000))
best_gensa <- result_gensa$par
best_gensa



