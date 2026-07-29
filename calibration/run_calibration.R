# calibration/run_calibration.R

# Searches for the best-fit values of param1, param2a, param2b via full
# factorial analysis.
#
# Output:
#   - calibration/calibration_results.RDS  (full grid results)
#   - Console output: best params to copy into calibrate_config.R

library(tidyverse)
library(parallel)


# ---- 1. SOURCE MODEL AND LOAD DATA ----------------------------------------

source("model/setup.R")               # defines l.inputs
source("model/helpers/source_all.R")  # loads all helper/module functions
source("model/simulation.R")          # loads f.run(), f.initialize()

# Load benchmark targets (benchmarking_helpers.R uses <<- to assign globals)
source("calibration/benchmarking_helpers.R")
bench_prev      <- benchmark_prev_by_age   # ages, condition, prev, ci_lo, ci_hi
bench_lifetable <- lifetable               # age, qx, rate

# Load microdata
microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")


# ---- 2. PARAMETER GRID -----------------------------------------------------

n_steps <- 6
grid <- expand.grid(
  param1  = seq(0.20, 2.85, length.out = n_steps),
  param2a = seq(1,    3,    length.out = n_steps),
  param2b = seq(1,    2.67, length.out = n_steps)
)
n_calib <- 10000
cat(sprintf("Grid: %d combinations | n = %d per run\n", nrow(grid), n_calib))

#


# ---- 3. GOF FUNCTION -------------------------------------------------------
# Computes normalized WSSD (per protocol sec. 2.3), over three equally weighted targets:
#   Total GOF = WSSD_prev_mci / n_prev_mci + WSSD_prev_dem / n_prev_dem + WSSD_mort / n_mort

compute_gof <- function(sim_output, bench_prev, bench_lifetable) {

  # -- Prevalence ---
  df_prev <- data.frame(
    age   = as.vector(sim_output[, "AGE",   ]),
    alive = as.vector(sim_output[, "ALIVE", ]),
    sev   = as.vector(sim_output[, "SEV",   ])
  ) %>%
    filter(alive == 1) %>%
    mutate(condition = case_when(
      sev == 0     ~ "mci",
      sev %in% 1:3 ~ "dem",
      TRUE         ~ "healthy"
    ))

  prev_by_age <- df_prev %>%
    filter(age %in% bench_prev$age) %>%
    group_by(age) %>%
    summarise(
      mci = sum(condition == "mci") / n(),
      dem = sum(condition == "dem") / n(),
      .groups = "drop"
    ) %>%
    pivot_longer(c(mci, dem), names_to = "condition", values_to = "model_prev")

  prev_compare <- bench_prev %>%
    mutate(SE = (ci_hi - ci_lo) / 3.92) %>%
    left_join(prev_by_age, by = c("age", "condition")) %>%
    filter(!is.na(model_prev), !is.na(SE), SE > 0)

  # Separate MCI and dementia so each is an equal-weight target (3 vs 6 points
  # would otherwise give dementia 2x the pull within a combined prevalence term)
  mci_compare <- prev_compare %>% filter(condition == "mci")
  dem_compare <- prev_compare %>% filter(condition == "dem")

  wssd_mci <- sum((mci_compare$model_prev - mci_compare$prev)^2 / mci_compare$SE^2)
  wssd_dem <- sum((dem_compare$model_prev - dem_compare$prev)^2 / dem_compare$SE^2)
  n_mci    <- nrow(mci_compare)
  n_dem    <- nrow(dem_compare)

  # -- Mortality ---
  # Both sides are annual conditional probabilities of death: P(dies during the year |
  # alive at its start). That is the life table's native unit (`qx`), the model's native
  # input (m.lifetable holds probabilities), and what the simulation actually draws, so no
  # scale conversion is needed on either side.
  #
  # Computed per individual, pairing each person's age at the start of a cycle with whether
  # they died during it, the same way compare_mortality() does, then joined to the benchmark
  # on age. The join must be on age rather than on cycle position: a death recorded at cycle
  # t is drawn using AGE at t-1 (see f.update_ALIVE), and cycle number equals age only for a
  # cohort that starts at exactly the life table's first age with no spread.
  n_cycle <- dim(sim_output)[1]
  alive   <- sim_output[, "ALIVE", ]

  model_mort <- data.frame(
    age       = as.vector(sim_output[-n_cycle, "AGE", ]),  # age at start of cycle
    was_alive = as.vector(alive[-n_cycle, ]) == 1,
    died      = as.vector(alive[-1, ]) == 0                # died during the cycle
  ) %>%
    filter(was_alive) %>%
    group_by(age) %>%
    summarise(n_at_risk = n(),
              n_deaths  = sum(died, na.rm = TRUE),
              .groups   = "drop") %>%
    mutate(model_prob = n_deaths / n_at_risk)

  # The life table's last age carries qx = 1 (death is absorbing where the table ends), an
  # artifact of the table rather than a model target, so it is dropped explicitly instead of
  # being swallowed by na.rm while still counting toward n_mort. The n_at_risk floor guards
  # against ages thin enough that model_prob is mostly noise, which a high-mortality
  # parameter set can produce at the oldest ages.
  mort_compare <- bench_lifetable %>%
    filter(qx < 1) %>%
    inner_join(model_mort, by = "age") %>%
    filter(n_at_risk >= 20) %>%
    arrange(age)

  # SE for life table (no published CI): proportional SE = 5% of qx, floored at 1e-4
  se_mort   <- pmax(mort_compare$qx * 0.05, 1e-4)
  wssd_mort <- sum((mort_compare$model_prob - mort_compare$qx)^2 / se_mort^2)
  n_mort    <- nrow(mort_compare)

  # -- Total: three equal-weight targets (MCI prev, dementia prev, mortality) --
  gof <- wssd_mci / n_mci + wssd_dem / n_dem + wssd_mort / n_mort

  list(gof = gof, wssd_mci = wssd_mci, wssd_dem = wssd_dem, wssd_mort = wssd_mort,
       n_mci = n_mci, n_dem = n_dem, n_mort = n_mort)
}


# ---- 4. WORKER FUNCTION (one grid point) -----------------------------------

run_one <- function(row, l.inputs, microdata, bench_prev, bench_lifetable, n_calib) {

  inputs_local <- l.inputs
  inputs_local[["n.ind"]]   <- n_calib
  inputs_local[["n.cycle"]] <- 51

  # Calibration-specific base settings (mirrors calibrate_config.R)
  inputs_local[["p.HCARE_start"]]   <- c(0.25, 0.75)
  inputs_local[["hr.mort_mci_age"]] <- c(1, 1, 1)
  inputs_local[["hr.mort_mod_age"]] <- c(1, 1, 1)
  inputs_local[["hr.mort_sev_age"]] <- c(1, 1, 1)
  inputs_local[["seed_stochastic"]] <- 20250624

  # Apply search parameters
  r.CDRslow_base <- l.inputs[["r.CDRslow_mean"]]  # base value (0.6) from setup.R
  inputs_local[["param1"]]         <- row[["param1"]]
  inputs_local[["param2a"]]        <- row[["param2a"]]
  inputs_local[["param2b"]]        <- row[["param2b"]]
  inputs_local[["m.hr_mci"]]       <- l.inputs[["m.hr_mci"]] * row[["param1"]]
  inputs_local[["r.CDRslow_mean"]] <- (seq(0, 1, length.out = 51)^row[["param2a"]]) *
                                        (row[["param2b"]] * r.CDRslow_base)

  # Run simulation (f.run directly; skip figure generation). f.out_aggregate() is not called:
  # every GOF target is computed from the raw a.out array, so aggregating would build 181 list
  # elements per grid point (~10% of each run's time) for nothing. Re-add it here if a future
  # target needs an aggregated quantity such as reside time or age at onset.
  result <- tryCatch({
    output     <- f.run(l.inputs = inputs_local, microdata = microdata, printLevel = 0)
    gof_vals   <- compute_gof(output, bench_prev, bench_lifetable)
    # Return a flat named list of scalars so bind_rows() works cleanly
    list(
      param1    = row[["param1"]],
      param2a   = row[["param2a"]],
      param2b   = row[["param2b"]],
      gof       = gof_vals$gof,
      wssd_mci  = gof_vals$wssd_mci,
      wssd_dem  = gof_vals$wssd_dem,
      wssd_mort = gof_vals$wssd_mort,
      n_mci     = gof_vals$n_mci,
      n_dem     = gof_vals$n_dem,
      n_mort    = gof_vals$n_mort,
      error     = FALSE,
      error_msg = NA_character_
    )
  }, error = function(e) {
    list(
      param1    = row[["param1"]],
      param2a   = row[["param2a"]],
      param2b   = row[["param2b"]],
      gof       = Inf,
      wssd_mci  = NA_real_,
      wssd_dem  = NA_real_,
      wssd_mort = NA_real_,
      n_mci     = NA_integer_,
      n_dem     = NA_integer_,
      n_mort    = NA_integer_,
      error     = TRUE,
      error_msg = conditionMessage(e)
    )
  })

  result
}


# ---- 5-7. CALIBRATION RUN (parallel execution, save, best-fit + boundary check) --
# Wrapped as a function so a second pass can be run with a different (e.g. finer)
# grid without duplicating this logic. Uses l.inputs, microdata, bench_prev,
# bench_lifetable, compute_gof, run_one from the enclosing script.

run_calibration_grid <- function(grid, n_calib, n_steps,
                                  results_file = "calibration/calibration_results.RDS") {

  # -- Parallel execution --
  n_cores <- max(1, detectCores() - 1)
  cat(sprintf("Starting cluster: %d cores\n", n_cores))
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)   # clean up even if function errors

  clusterSetRNGStream(cl, iseed = 20250624)

  # Load packages on workers
  clusterEvalQ(cl, {
    library(tidyverse)
    library(rlang)
    library(scales)
    library(patchwork)
    library(flextable)
  })

  # Source function files on workers (explicit list avoids benchmarking_helpers.R side effects)
  clusterEvalQ(cl, {
    source("model/helpers/run_wrappers.R")
    source("model/helpers/output_formatters.R")
    source("model/helpers/epi_helpers.R")
    source("model/helpers/generic_helpers.R")
    source("model/helpers/plotting_helpers.R")
    source("model/helpers/program_utility_helpers.R")
    lapply(list.files("model/modules", pattern = "^module_.*\\.R$", full.names = TRUE), source)
    lapply(list.files("model/config",  pattern = "_config\\.R$",    full.names = TRUE), source)
    source("model/simulation.R")
  })

  # Export data and functions needed by workers
  clusterExport(cl, c("l.inputs", "microdata", "bench_prev", "bench_lifetable",
                      "grid", "n_calib", "compute_gof", "run_one"))

  cat(sprintf("Running %d combinations on %d cores...\n", nrow(grid), n_cores))
  t_start <- proc.time()

  results_list <- parLapply(cl, seq_len(nrow(grid)), function(i) {
    run_one(as.list(grid[i, ]), l.inputs, microdata, bench_prev, bench_lifetable, n_calib)
  })

  t_elapsed <- proc.time() - t_start
  cat(sprintf("Done. Elapsed: %.1f minutes\n", t_elapsed["elapsed"] / 60))

  # -- Collect and save results --
  results <- bind_rows(lapply(results_list, as_tibble))
  saveRDS(results, results_file)
  cat(sprintf("Results saved to %s\n", results_file))

  n_errors <- sum(results$error, na.rm = TRUE)
  if (n_errors > 0) {
    cat(sprintf("WARNING: %d combinations failed (see results$error_msg)\n", n_errors))
  }

  # -- Best params and boundary check --
  best <- results %>%
    filter(!error) %>%
    slice(which.min(gof))

  cat("\n=== BEST PARAMETER SET ===\n")
  cat(sprintf("  param1  = %.4f  (searched %.2f – %.2f)\n",
              best$param1,  min(grid$param1),  max(grid$param1)))
  cat(sprintf("  param2a = %.4f  (searched %.2f – %.2f)\n",
              best$param2a, min(grid$param2a), max(grid$param2a)))
  cat(sprintf("  param2b = %.4f  (searched %.2f – %.2f)\n",
              best$param2b, min(grid$param2b), max(grid$param2b)))
  cat(sprintf("  GOF     = %.4f  (WSSD_mci/n = %.4f | WSSD_dem/n = %.4f | WSSD_mort/n = %.4f)\n",
              best$gof,
              best$wssd_mci  / best$n_mci,
              best$wssd_dem  / best$n_dem,
              best$wssd_mort / best$n_mort))

  # Boundary check
  step_sizes <- c(
    param1  = diff(range(grid$param1))  / (n_steps - 1),
    param2a = diff(range(grid$param2a)) / (n_steps - 1),
    param2b = diff(range(grid$param2b)) / (n_steps - 1)
  )
  at_boundary <- c(
    param1  = (best$param1  - min(grid$param1)  < step_sizes["param1"]  * 0.01) |
              (max(grid$param1)  - best$param1  < step_sizes["param1"]  * 0.01),
    param2a = (best$param2a - min(grid$param2a) < step_sizes["param2a"] * 0.01) |
              (max(grid$param2a) - best$param2a < step_sizes["param2a"] * 0.01),
    param2b = (best$param2b - min(grid$param2b) < step_sizes["param2b"] * 0.01) |
              (max(grid$param2b) - best$param2b < step_sizes["param2b"] * 0.01)
  )

  if (any(at_boundary)) {
    cat(sprintf("\n*** WARNING: boundary solution for: %s ***\n",
                paste(names(at_boundary)[at_boundary], collapse = ", ")))
    cat("Expand the range for that parameter and re-run.\n")
  } else {
    cat("\nInterior optimum — no boundary issues.\n")
  }

  list(results = results, best = best)
}

# ---- Run pass 1 (coarse grid) ----------------------------------------------

calib_1  <- run_calibration_grid(grid, n_calib, n_steps, 
                                 results_file = "calibration/calibration_results_20260728.RDS")
results  <- calib_1$results
best     <- calib_1$best

# ---- Run pass 2 (finer grid) -----------------------------------------------
n_steps <- 10
grid <- expand.grid(
  param1  = seq(best$param1*(1-0.33), best$param1*(1+0.33), length.out = n_steps),
  param2a = seq(best$param2a*(1-0.33), best$param2a*(1+0.33), length.out = n_steps),
  param2b = seq(best$param2b*(1-0.33), best$param2b*(1+0.33), length.out = n_steps)
)
n_calib <- 10000

calib_2 <- run_calibration_grid(grid, n_calib, n_steps,
                                results_file = "calibration/calibration_results_20260728_pass2.RDS")
results <- calib_2$results
best    <- calib_2$best


# ---- 8. COPY-PASTE OUTPUT FOR calibrate_config.R --------------------------

cat("\n=== UPDATE calibrate_config.R WITH: ===\n")
cat(sprintf('  inputs[["param1"]]  <- %.4f\n', best$param1))
cat(sprintf('  inputs[["param2a"]] <- %.4f\n', best$param2a))
cat(sprintf('  inputs[["param2b"]] <- %.4f\n', best$param2b))


# ---- 9. GOF SURFACE PLOTS --------------------------------------------------
# For each parameter, fix the other two at their best values and plot GOF.
# A smooth bowl shape confirms no coarseness issues; a flat edge suggests
# the grid may need refinement in that direction.

plot_gof_surface <- function(results, best, free_param, fixed_params) {
  results %>%
    filter(
      .data[[fixed_params[1]]] == best[[fixed_params[1]]],
      .data[[fixed_params[2]]] == best[[fixed_params[2]]]
    ) %>%
    ggplot(aes(x = .data[[free_param]], y = gof)) +
    geom_line(color = "steelblue") +
    geom_point(color = "steelblue", size = 2) +
    geom_point(data = best, aes(x = .data[[free_param]], y = gof),
               color = "red", size = 3) +
    labs(
      title  = sprintf("GOF surface: %s", free_param),
      subtitle = sprintf("%s = %.4f, %s = %.4f (fixed at best)",
                         fixed_params[1], best[[fixed_params[1]]],
                         fixed_params[2], best[[fixed_params[2]]]),
      x = free_param, y = "GOF"
    ) +
    theme_minimal(base_size = 13)
}

fig_p1  <- plot_gof_surface(results, best, "param1",  c("param2a", "param2b"))
fig_p2a <- plot_gof_surface(results, best, "param2a", c("param1",  "param2b"))
fig_p2b <- plot_gof_surface(results, best, "param2b", c("param1",  "param2a"))

fig_gof_surface <- fig_p1 / fig_p2a / fig_p2b +
  plot_annotation(title = "GOF surface by calibration parameter",
                  subtitle = "Red dot = best-fit point | Smooth bowl = grid is adequate")
print(fig_gof_surface)
