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

# setup.R holds UNCALIBRATED defaults. calibrate() overlays the calibrated
# values from model/config/calibrate_config.R and sets n.ind / n.cycle.
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 10000)

microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
results <- f.wrap_run(l.inputs_calibrated, microdata = microdata)
```

There is no automated test suite. Calibration is run via `calibration/run_calibration.R`, and validation via analysis-specific scripts.

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

- `model/setup.R` — all model parameters (`l.inputs`), uncalibrated defaults
- `model/config/calibrate_config.R` — `calibrate()` and the calibrated parameter
  values; the single source of truth for these, never restate them in docs
- `model/simulation.R` — main simulation engine (`f.run`, `f.initialize`)
- `model/modules/MODULES.md` — module-level documentation
- `calibration/run_calibration.R` — full factorial calibration grid search
- `calibration/GRAM_calibration_validation_protocol.docx` — calibration and
  validation protocol (the working copy; there is no Markdown twin)
- `analyses/testing_strategies/` — active analysis (Paper 2: Testing Strategies)
- `docs/ONBOARDING.md` — guide for country-specific adaptations
- `archive/` — superseded files, kept for the record; nothing here is live

## Git Workflow

- `main`: stable releases, `develop`: active development
- Country adaptations: `country/<name>`, features: `feature/<name>`, fixes: `fix/<name>`
- Core model code in `model/` should not be modified for country adaptations; use config overrides instead
