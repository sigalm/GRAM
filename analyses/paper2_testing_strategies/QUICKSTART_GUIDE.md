# Quickstart Guide: Running Alternative Testing Strategy Scenarios

This guide will help you run simulations and generate results for different cognitive testing strategies in the GRAM model.

---

## Overview

The testing strategies analysis evaluates different approaches to cognitive screening:
- **Inclusive (Universal)**: Regular testing regardless of cognitive concerns
- **Selective**: Testing based on prompted cognitive concerns
- **Reactive**: Testing only when concerns are brought up spontaneously

Each strategy can be configured with different testing frequencies and follow-up protocols.

---

## File Structure

```
analyses/paper2_testing_strategies/
├── gram_paper2_testingstrategies.R     # Main analysis script
├── gram_paper2_testing_pcp_fu.R        # PCP follow-up analysis
├── gram_evolution_charts.R             # Helper functions for early diagnosis analysis
├── test_performance_helpers.R          # Functions for test performance metrics & plotting
├── bha_scenarios/                      # Scenario configuration files
│   ├── scenario_TEMPLATE_config.R          # Template for creating new scenarios
│   ├── scenario_u3bhapos_rand50_config.R   # Universal every 3y, 50% random uptake
│   ├── scenario_u3bhapos_nonrand_config.R  # Universal every 3y, non-random uptake
│   ├── scenario_s1bhapos_emr_config.R      # Selective annually, EMR-triggered
│   ├── scenario_s1bhapos_question_config.R # Selective annually, question-triggered
│   ├── scenario_r1bhapos_config.R          # Reactive annually
│   ├── scenario_u1hybrid_config.R          # Hybrid strategy
│   └── scenario_*pcp*_config.R             # Variants with PCP follow-up
├── planning/                           # Planning documents
├── validation/                         # Validation outputs
├── sim_results/                        # Simulation output files (.rds)
├── test_perf_results/                  # Processed test performance data
└── plots/                              # Generated figures
```

Run `list.files("analyses/paper2_testing_strategies/bha_scenarios")` for the current set —
scenarios are added and retired often enough that this listing goes out of date.

---

## Quick Start: Running Existing Scenarios

### Step 1: Open the Main Script

Open `gram_paper2_testingstrategies.R` in RStudio.

### Step 2: Run Setup and Calibration

```r
# Load required libraries and scripts
source("model/setup.R")
source("model/simulation.R")
source("model/helpers/source_all.R")   # also loads calibrate() from model/config/
library(tableone)

# Load microdata
sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")

# Apply the calibrated parameter values and set the cohort size
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 100000)
```

`calibrate()` is defined in `model/config/calibrate_config.R` and loaded by
`source_all.R` — it does not come from `calibration/benchmarking_helpers.R`, which only
supplies benchmark targets and comparison functions. It applies the calibrated `param1`,
`param2a` and `param2b`, and sets `n.ind` to `n` and `n.cycle` to 51; it does not re-run
the calibration search (that is `calibration/run_calibration.R`), so it returns quickly.

### Step 3: Select Scenarios to Run

The script defines scenarios in a list. Modify which scenarios to run by adjusting the indices:

```r
# The three main strategies (u3bhapos_rand50, s1bhapos_emr, r1bhapos)
for (scen in names(scenario_list[1:3])) {
  # ... simulation code ...
}

# All scenarios in the list, including the sensitivity analyses
for (scen in names(scenario_list)) {
  # ... simulation code ...
}

# To run specific scenarios only:
for (scen in c("u3bhapos", "r1bhapos")) {
  # ... simulation code ...
}
```

### Step 4: Run Simulations

Execute the simulation loop. Each scenario will:
- Load its configuration file
- Run the simulation
- Save results with a timestamp to `sim_results/`

**Note**: Each simulation can take several minutes to hours depending on sample size.

### Step 5: Generate Test Performance Plots

After simulations complete, run the evolution charts section:

```r
# Process simulation outputs
for (scen in names(scenario_list[1:3])) {
  output <- latest_rds(scen)$output
  test_data <- post_processing_outputs(output)
  saveRDS(test_data, file = file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, ".rds")))
}

# Combine and plot results
test_data_combined <- data.frame()
for (i in seq_along(names(scenario_list[1:3]))) {
  scen <- names(scenario_list)[i]
  test_data <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", paste0(scen, ".rds"))) %>%
    mutate(scenario = scen)
  test_data_combined <- rbind(test_data_combined, test_data)
}

# Create plots
subtitles <- c("Inclusive testing, every 3 years",
               "Selective testing, annual",
               "Reactive testing, annual")
names(subtitles) <- names(scenario_list)[1:3]

plot_test_results(test_data_combined, ages = 65:80, show_early_pos = FALSE, 
                  scenario_names = subtitles, y_max = 75000)

plot_testers(test_data_combined, ages = 65:80, scenario_names = subtitles)

ggsave("analyses/paper2_testing_strategies/plots/no-early-positives.jpeg", height = 10, width = 8)
```

---

## Creating a New Scenario

### Step 1: Copy the Template

1. Navigate to `bha_scenarios/`
2. Copy `scenario_TEMPLATE_config.R`
3. Rename to `scenario_<YOURNAME>_config.R`

### Step 2: Configure Your Scenario

Edit the new file to define your testing strategy:

```r
scenario_inputs <- list(
  
  # ============================================================================
  # SCENARIO IDENTIFICATION
  # ============================================================================
  title       = as.factor("YOURNAME"),
  description = "Brief description of your testing scenario",

  # ============================================================================
  # TEST PARAMETERS
  # ============================================================================
  test        = "BHA-GS",
  sensitivity = l.inputs[["sens_BHAGS"]],
  specificity = l.inputs[["spec_BHAGS"]],
  
  # ============================================================================
  # GLOBAL PARAMETERS
  # ============================================================================
  HCARE = 1,    # 1 = requires healthcare access, 0 = ignore
  DX    = 0,    # 0 = exclude those with prior diagnosis, 1 = include all
  
  # ============================================================================
  # CORE SCENARIO PARAMETERS
  # ============================================================================
  age_first_test  = 65,     # Starting age for testing
  age_stop_test   = 80,     # Age to stop testing
  
  # Cognitive concern probabilities - choose one:
  probs_cogcon = l.inputs[["m.cogcon"]],              # Universal (everyone tested)
  # probs_cogcon = l.inputs[["m.cogcon_selective"]],  # Selective (prompted concerns)
  # probs_cogcon = l.inputs[["m.cogcon_reactive"]],   # Reactive (spontaneous concerns)
  
  prob_pcpfu      = NULL,   # Probability of PCP follow-up (NULL = no follow-up)
  repeat_interval = 1,      # Years between tests
  cohort_split    = NULL,   # Split cohort (e.g., 3 = test 1/3 each year)
  
  # ============================================================================
  # STOP RULE
  # ============================================================================
  stop_rule = function(...) {
    args <- list(...)
    args$any_BHA_pos == TRUE  # Stop after first positive
  }
)
```

### Step 3: Add to Scenario List

In `gram_paper2_testingstrategies.R`, add your scenario to the list:

```r
scenario_list <- list(
  u3bhapos_rand50   = "scenario_u3bhapos_rand50_config.R",
  s1bhapos_emr      = "scenario_s1bhapos_emr_config.R",
  r1bhapos          = "scenario_r1bhapos_config.R",
  s1bhapos_question = "scenario_s1bhapos_question_config.R",
  u3bhapos_nonrand  = "scenario_u3bhapos_nonrand_config.R",
  yourname          = "scenario_yourname_config.R"  # Add your scenario here
)
```

### Step 4: Run Your Scenario

Update the loop to include your scenario (or provide index of scenario to run):

```r
for (scen in names(scenario_list)) {
  # ... simulation code runs all scenarios including yours ...
}
```

---

## Key Configuration Options

### Testing Frequency

- **`repeat_interval`**: Years between tests (1 = annual, 3 = every 3 years)
- **`cohort_split`**: Divide cohort into groups tested in rotation
  - Example: `cohort_split = 3` with `repeat_interval = 3` means 1/3 of cohort tested each year, with each group tested once every 3 years

### Testing Eligibility

- **`probs_cogcon`**: Determines who gets tested based on cognitive concerns
  - `m.cogcon`: Universal testing (everyone eligible)
  - `m.cogcon_selective`: Selective testing (prompted concerns)
  - `m.cogcon_reactive`: Reactive testing (spontaneous concerns only)
The base model selection criteria was based on whether the individual had a cognitive concern (m.cogcon > 0). Functionally, this variable is a vector of probabilities for getting selected for testing based on true cognitive status. Therefore, it can be used to implement any probabilistic selection criteria, simply replace the values of m.cogcon with your own selection probabilities.

### Age Range

- **`age_first_test`**: Age at which individuals are eligibile for testing
- **`age_stop_test`**: Age at which testing ceases

### Stop Rules

Define when to stop testing an individual. The stop_rule function is called within `f.update_BHA()` in the medical record module (`model/modules/module_medical_record.R`).

**Available Arguments:**

The stop_rule function receives the following arguments (all are vectors with one value per individual):

- **`any_BHA_pos`**: Logical indicating if the individual ever had a positive BHA result
- **`NP`**: Neuropsychological test result from previous cycle (lagged)
- **`PET`**: PET imaging result from previous cycle (lagged)
- **`any_PCP_pos`**: Logical indicating if the individual ever had a positive PCP evaluation
- **`repeat_after_FP`**: Scenario parameter for pausing testing after false positives

**Examples:**

```r
stop_rule = function(...) {
  args <- list(...)
  
  # Stop after first positive BHA
  args$any_BHA_pos == TRUE
  
  # Or continue testing regardless (never stop)
  # FALSE
  
  # Or stop after PCP confirmation
  # args$any_PCP_pos == TRUE
  
  # Or stop after neuropsych testing
  # args$NP == 1
}
```

**Adding Custom Arguments:**

If you need additional arguments for your stop rule, you must update `f.update_BHA()` in `model/modules/module_medical_record.R` (around line 162-166) to pass the new arguments to the stop_rule function call. 

---

## Understanding Results

### Simulation Output

Each simulation saves an `.rds` file containing a 3D array:
- **Dimension 1**: Time cycles
- **Dimension 2**: Variables/attributes (cognitive states, test results, diagnoses)
- **Dimension 3**: Individuals

### Test Performance Metrics

The `post_processing_outputs()` function calculates:
- **TP (True Positives)**: Correctly identified impaired individuals
- **FP (False Positives)**: Healthy individuals incorrectly flagged
- **TN (True Negatives)**: Correctly identified healthy individuals
- **FN (False Negatives)**: Missed impaired individuals
- **Early Positives**: Healthy individuals who later develop impairment
- **Non-testers**: Eligible individuals who were never tested due to selection probabilities

### Visualization

Two main plot types are generated:

1. **`plot_test_results()`**: Shows cumulative counts of test outcomes by cognitive status
   - Red lines = Impaired individuals
   - Blue lines = Healthy individuals
   - Filled circles = Positive test results
   - Empty circles = Negative test results
   - X marks = Not tested

2. **`plot_testers()`**: Shows testing coverage over time
   - Tested vs. not tested by cognitive status

---

## Helper Functions

### Finding Latest Results

```r
# Get most recent simulation output for a scenario
latest_result <- latest_rds("u3bhapos")
output_array <- latest_result$output
```

### Analyzing Specific Ages

```r
# Extract data for age 65
age_65_data <- as.data.frame(t(output_array[65-50+1,,]))

# Calculate prevalence among undiagnosed at age 65
prev_in_undx_65 <- age_65_data %>%
  filter(ALIVE == 1, DX == 0) %>%
  mutate(status = case_when(
    SYN < 1 ~ "healthy",
    SEV == 0 ~ "mci",
    SEV >= 1 ~ "dementia"
  )) %>%
  group_by(status) %>%
  summarise(n = n()) %>%
  mutate(prev = n/sum(n))
```

### Custom Plotting

```r
# Plot cumulative diagnoses
plot_cumulative_count(
  output_array,
  variables = list(
    "SYN" = list(variable_name = "SYN", condition_value = 1),
    "DX" = list(variable_name = "DX", condition_value = 1)
  ),
  plot_title = "Cumulative Diagnoses",
  scenario_name = "My Scenario"
)
```

---

## Troubleshooting

### Simulation Takes Too Long

- Reduce sample size in calibration: `calibrate(inputs = l.inputs, n = 50000)`
- Run fewer scenarios at once
- Use a smaller age range in plotting functions
- Avoid re-running the same scenario to save time. You can load the results from the `sim_results/` directory then use `post_processing_outputs()` to calculate test performance metrics. For large simulations (e.g., 100,000 individuals), post-processing can also take a long time. Post-processed outputs can be saved to the `test_perf_results/` directory and re-loaded for plotting and analysis.

### Memory Issues

The script includes garbage collection between scenarios:
```r
rm(config, result)
invisible(gc())
```

If memory issues persist, run scenarios individually rather than in a loop, or consider using the memory limit of your R session.

### Results Not Found

Check that:
1. Simulations completed successfully (check `sim_results/` for `.rds` files)
2. Scenario names match exactly between `scenario_list` and file names
3. File paths are correct relative to project root

---

## Example Workflow

Here's a complete workflow for comparing two testing strategies:

```r
# 1. Setup -- source_all.R is required: it defines calibrate() and load_scenario()
source("model/setup.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")
source("model/helpers/source_all.R")
source("analyses/paper2_testing_strategies/test_performance_helpers.R")
library(tableone)
library(flextable)
sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")

# 2. Apply calibrated parameters
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 100000)  # you only need to run this once

# 3. Define scenarios
scenario_list <- list(
  u3bhapos_rand50 = "scenario_u3bhapos_rand50_config.R",
  s1bhapos_emr    = "scenario_s1bhapos_emr_config.R"
)

# 4. Run simulations
for (scen in names(scenario_list)) {
  config_file <- file.path("analyses/paper2_testing_strategies/bha_scenarios", 
                           scenario_list[[scen]])
  config <- load_scenario(config_file, l.inputs_calibrated)
  result <- f.wrap_run(config, microdata = sample1)
  datetime_suffix <- format(Sys.time(), "%Y%m%d_%H%M%S")
  saveRDS(result, file = file.path("analyses/paper2_testing_strategies/sim_results", 
                                   paste0("scenario_", scen, "_sim_", datetime_suffix, ".rds")))
  rm(config, result)
  invisible(gc())
}

# 5. Process results
for (scen in names(scenario_list)) {
  output <- latest_rds(scen)$output
  test_data <- post_processing_outputs(output)
  saveRDS(test_data, file = file.path("analyses/paper2_testing_strategies/test_perf_results", 
                                      paste0(scen, ".rds")))
}

# 6. Generate plots
test_data_combined <- data.frame()
for (scen in names(scenario_list)) {
  test_data <- readRDS(file.path("analyses/paper2_testing_strategies/test_perf_results", 
                                 paste0(scen, ".rds"))) %>%
    mutate(scenario = scen)
  test_data_combined <- rbind(test_data_combined, test_data)
}

subtitles <- c("Universal testing every 3 years", "Selective testing annually")
names(subtitles) <- names(scenario_list)

plot_test_results(test_data_combined, ages = 65:80, show_early_pos = FALSE,
                  scenario_names = subtitles, y_max = 75000)
ggsave("analyses/paper2_testing_strategies/plots/my_comparison.jpeg", 
       height = 10, width = 8)
```

---

## Additional Resources

- **Model parameters**: See `docs/methods_appendix.Rmd` for detailed parameter descriptions
- **Calibration benchmarks**: See `calibration/GRAM calibration benchmarks.xlsx`
- **Cognitive concern data**: See `data/cogcon/` for probability matrices
- **Planning documents**: See `planning/` folder for scenario design rationale

---

## Questions?

If you encounter issues or need to modify the analysis in ways not covered here, review:
1. The scenario template for all available configuration options
2. The test performance helper functions for custom metrics
3. The main simulation script and medical record module for the overall workflow structure
