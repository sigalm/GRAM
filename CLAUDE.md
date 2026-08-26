# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GRAM (GBHI Resource Allocation Model) is an R-based microsimulation model that simulates cognitive impairment progression and intervention evaluation in aging populations. It models individual trajectories through cognitive health states: Healthy → MCI → Mild/Moderate/Severe Dementia → Death.

## Running the Model

```r
# Open GRAM.Rproj in RStudio (sets working directory), then:
source("model/setup.R")              # Load l.inputs (all parameters)
source("model/helpers/source_all.R") # Load helpers, modules & configs
source("model/simulation.R")         # Load simulation functions

# setup.R already applies the calibrated parameters, read from
# model/config/calibrated_params.R. There is nothing to call: l.inputs is
# calibrated as sourced. Set n.ind / n.cycle for the analysis at hand.
l.inputs[["n.ind"]] <- 10000

microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
results <- f.wrap_run(l.inputs, microdata = microdata)
```

There is no automated test suite. Calibration is run via `calibration/run_calibration.R`, and validation via analysis-specific scripts.

## Before Calling Something a Bug

Outcome definitions in this repo encode deliberate analytic intent that is not
recoverable from the code alone. Two series that look inconsistent with each
other are usually answering different questions on purpose.

When you notice what looks like a bug, an inconsistency, or a definition that
"should" be changed: **ask what the quantity is meant to measure before
asserting anything is wrong.** State the observation, ask the question, and wait
for the answer before proposing a fix. One question up front is far cheaper than
an exchange spent working backwards from a wrong assumption about the goal.

Be especially slow to assert on:

- numerator and denominator choices in performance metrics
- whether a series is a stock (who is in this state at time t) or a flow
  (how many events happened in cycle t)
- exit rules: when someone stops being counted, and why
- an eligibility gate that appears on one series but not on a sibling

Worked example, so this one is not re-litigated. In
`analyses/testing_strategies/test_performance_helpers.R`, `fn` counts a person
in every cycle after a negative test and carries no eligibility gate, while its
sibling `notest_fn` gates on lagged DX. That asymmetry is intentional. `fn` is a
**verdict**: the test looked at an impaired person and said no, and that error
stands until a later test overturns it — losing the chance to be corrected (by a
clinical DX ending eligibility) does not unmake it. `notest_fn` is a **coverage
gap**: an eligible impaired person the program has not reached, which ceases to
exist once they are no longer eligible, because there was never a verdict to be
wrong about. A verdict persists; a gap is transient. The intended quantity is
errors outstanding among the living: corrections leave via `cummax`, deaths
leave because `SYN` is `NA` after death, and DX carryover stays.

## Architecture

### Three-Dimensional Array Core

The simulation tracks state in 3D arrays: `a.out[cycle, attribute, individual]` and `a.random[cycle, attribute, individual]` (pre-generated random numbers for deterministic reproducibility).

### Simulation Loop (`f.run()` in `model/simulation.R`)

Each cycle (year), modules execute in order:
1. `f.module_mortality()` — update ALIVE status
2. `f.module_socdem()` — update AGE, HCARE
3. `f.module_true_health()` — update SYN, SEV, CDR (true cognitive status)
4. `f.module_medical_record()` — update BHA, CDR_obs, DX (observed/tested)
5. `f.module_treatment()` — update TX, TX2, LTC

`f.wrap_run()` (in `model/helpers/run_wrappers.R`) orchestrates: run simulation → aggregate results → generate figures.

### Module Pattern

Each module in `model/modules/` follows the same structure:
- Module wrapper `f.module_<name>()` orchestrates updates for a cycle
- Update functions `f.update_<attribute>()` modify specific attributes
- All use lagged values (`.lag`) from t-1 and write to current cycle t
- Only living individuals are updated (`a.out[t,"ALIVE",]==1`)

### Scenario Configuration

Scenarios are R lists in `analyses/*/bha_scenarios/*_config.R` that override base `l.inputs` parameters. They define test type, sensitivity/specificity, age ranges, frequency, eligibility, and stop rules.

## Naming Conventions

| Prefix | Meaning | Example |
|--------|---------|---------|
| `f.` | Function | `f.update_SYN()` |
| `v.` | Vector | `v.AGE.lag` |
| `m.` | Matrix | `m.lifetable` |
| `l.` | List | `l.inputs` |
| `a.` | Array | `a.out` |
| `n.` | Number / count | `n.cycle` |
| `p.` | Probability | `p.APOE4_start` |
| `r.` | Rate | `r.CDRslow_mean` |
| `hr.` / `rr.` | Hazard ratio / relative risk | `hr.mort_mci` |
| `u.` / `c.` | Utility / cost | `u.mci`, `c.bha` |
| `log_` | Natural log of a coefficient | `log_APOE4` |
| `.lag` | Lagged value (t-1) | `v.SYN.lag` |

Attributes are mostly integer-encoded — ALIVE 0/1, SEX 1=male/2=female,
SEV 0-3, BHA -9=not tested/0=neg/1=pos, DX/APOE4/HCARE/TX/TX2/LTC 0/1.
Two exceptions worth knowing:

- `SYN` takes 0 / **0.5** / 1, where 0.5 is transitional cognitive impairment
  (TCI), a two-cycle tunnel state tracked alongside it by the `TCI` attribute.
- `EDU` is **years** of education, not an ordinal level.

`v.attr_names` in `model/setup.R` is the authoritative list; it currently holds
36 attributes, including cost accumulators (`COST_test`, `COST_fu`, `COST_tx2`,
`COST_care`, `COST_tx`) and test results (`PCP`, `PET`, `NP`).

## Key Files

- `model/setup.R` — all model parameters (`l.inputs`); sources the calibrated
  values at the end, so `l.inputs` is calibrated as sourced
- `model/config/calibrated_params.R` — GENERATED. The calibrated parameter values
  and their provenance; the single source of truth for these, never restate them
  in docs or analysis scripts. Regenerate it, don't edit it
- `model/simulation.R` — main simulation engine (`f.run`, `f.initialize`)
- `model/modules/MODULES.md` — module-level documentation
- `calibration/run_calibration.R` — full factorial calibration grid search;
  rewrites `calibrated_params.R` with the pass-2 best fit when it finishes
- `calibration/write_calibrated_params.R` — `f.write_calibrated_params()`, which
  promotes a saved results file into `calibrated_params.R`
- `docs/01_natural_history_supplement.Rmd` — model structure, parameters and
  their sources
- `docs/02_calibration_validation_supplement.Rmd` — calibration and validation
  protocol: targets, goodness-of-fit measure, acceptance criteria
- `analyses/testing_strategies/` — active analysis (Paper 2: Testing Strategies)
- `docs/ONBOARDING.md` — guide for country-specific adaptations
- `archive/` — superseded files, kept for the record; nothing here is live

## Git Workflow

- `main`: stable releases, `develop`: active development
- Country adaptations: `country/<name>`, features: `feature/<name>`, fixes: `fix/<name>`
- Core model code in `model/` should not be modified for country adaptations; use config overrides instead
