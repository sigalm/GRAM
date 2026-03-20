# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

GRAM (GBHI Resource Allocation Model) is an R-based microsimulation model that simulates cognitive impairment progression and intervention evaluation in aging populations. It models individual trajectories through cognitive health states: Healthy → MCI → Mild/Moderate/Severe Dementia → Death.

## Running the Model

```r
# Open GRAM.Rproj in RStudio (sets working directory), then:
source("model/setup.R")              # Load l.inputs (all parameters)
source("model/helpers/source_all.R") # Load helpers & modules
source("model/simulation.R")         # Load simulation functions

microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
results <- f.wrap_run(l.inputs, microdata = microdata)
```

There is no automated test suite. Validation is done via calibration notebooks (`calibration/benchmarking_v3.Rmd`) and analysis-specific validation scripts.

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

`f.wrap_run()` orchestrates: run simulation → aggregate results → generate figures.

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
| `.lag` | Lagged value (t-1) | `v.SYN.lag` |

Attributes are integer-encoded (e.g., ALIVE: 0/1, SEX: 1=male/2=female, SYN: 0=healthy/1=impaired, SEV: 0-3, BHA: -9=not tested/0=neg/1=pos).

## Key Files

- `model/setup.R` — all model parameters (`l.inputs`)
- `model/simulation.R` — main simulation engine (`f.run`, `f.initialize`, `f.wrap_run`)
- `model/modules/MODULES.md` — module-level documentation
- `analyses/paper2_testing_strategies/` — active analysis (Paper 2: Testing Strategies)
- `docs/ONBOARDING.md` — guide for country-specific adaptations

## Git Workflow

- `main`: stable releases, `develop`: active development
- Country adaptations: `country/<name>`, features: `feature/<name>`, fixes: `fix/<name>`
- Core model code in `model/` should not be modified for country adaptations; use config overrides instead
