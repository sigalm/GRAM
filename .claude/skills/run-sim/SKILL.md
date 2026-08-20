---
name: run-sim
description: Run a GRAM simulation with optional scenario configuration
argument-hint: "[scenario-config-file]"
disable-model-invocation: true
allowed-tools: Bash, Read, Glob
---

# Run GRAM Simulation

Run the GRAM microsimulation model. If a scenario config file is provided as an argument, run that scenario. Otherwise, run the base model with default parameters.

## Steps

1. **Source the model** in R:
   ```r
   source("model/setup.R")
   source("model/helpers/source_all.R")
   source("model/simulation.R")
   ```

2. **Load microdata**:
   ```r
   sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
   ```

3. **Run the simulation**:

   `l.inputs` is already calibrated when `setup.R` is sourced in step 1 — the values come
   from `model/config/calibrated_params.R`. There is no `calibrate()` call; only the cohort
   size needs setting. `load_scenario()` comes from `source_all.R` in step 1.

   - **If a scenario config file is provided** (`$ARGUMENTS`):
     Source the test performance helpers, set the cohort size, then load and run the scenario:
     ```r
     source("analyses/testing_strategies/test_performance_helpers.R")
     l.inputs[["n.ind"]] <- 100000
     config <- load_scenario("$ARGUMENTS", l.inputs)
     result <- f.wrap_run(config, microdata = sample1)
     ```
     Save output to `analyses/testing_strategies/sim_results/` with a datetime suffix.

   - **If no argument is provided**:
     Run the base model:
     ```r
     l.inputs[["n.ind"]] <- 100000
     result <- f.wrap_run(l.inputs, microdata = sample1)
     ```

4. **Execute** the above R code using `Rscript -e` via the Bash tool from the repository root directory. Combine all R commands into a single `Rscript -e` call. Use `cat()` to print key summary output so the user can see results.

5. **Report results** to the user: number of individuals simulated, number of cycles, and any summary statistics from `result$aggregated_results_totpop`.

## Notes

- The working directory must be the GRAM project root (where `GRAM.Rproj` lives).
- Scenario config files are typically located in `analyses/testing_strategies/bha_scenarios/`.
- Simulations with 100,000 individuals can take several minutes.
- If the user provides just a scenario name (e.g., `r1bhapos`, `u3bhapos_rand50`, `s1bhapos_emr`), look for the matching config file in `analyses/testing_strategies/bha_scenarios/` with the pattern `scenario_<name>_config.R`. Scenario names carry suffixes such as `_rand50`, `_nonrand`, `_emr` and `_question`, so list the directory rather than assuming a bare name resolves.
