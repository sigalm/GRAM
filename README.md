# GRAM: GBHI Resource Allocation Model

A microsimulation model for cognitive impairment progression and intervention evaluation in aging populations.

## Quick Start

```r
# 1. Open GRAM.Rproj in RStudio (sets working directory automatically)
#    Or manually: setwd("<path-to-GRAM>")

# 2. Source core scripts (in order)
source("model/setup.R")              # Loads l.inputs with parameters
source("model/helpers/source_all.R") # Loads helpers, modules, and calibrate()
source("model/simulation.R")         # Loads simulation functions

# 3. Apply the calibrated parameter values
#    setup.R holds uncalibrated defaults; calibrate() overlays the values in
#    model/config/calibrate_config.R and sets the cohort size and horizon.
l.inputs_calibrated <- calibrate(inputs = l.inputs, n = 10000)

# 4. Run a simulation
microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
results <- f.wrap_run(l.inputs_calibrated, microdata = microdata)
```

Skipping step 3 runs the model on uncalibrated defaults, which will not reproduce
published results.

## Directory Structure

```
GRAM/
├── model/                    # Core model code
│   ├── setup.R               # Model initialization and parameters (l.inputs)
│   ├── simulation.R          # Main simulation engine
│   ├── modules/              # Core simulation functions
│   │   └── MODULES.md        # Module documentation
│   ├── helpers/              # Utility functions (e.g., formatting, plotting)
│   │   └── source_all.R      # Sources all modules, helpers, and configs
│   └── config/
│       └── calibrate_config.R  # calibrate(): the calibrated parameter values
│
├── data/                     # Input data files
│   ├── acs_data/             # Population microdata (ACS-based)
│   ├── mortality/            # Life tables and mortality rates
│   ├── mci_incidence/        # MCI incidence rates by age
│   ├── cogcon/               # Cognitive concern probabilities
│   └── README.md             # Data provenance documentation
│
├── calibration/              # Calibration protocol, grid search, benchmarks
│   ├── run_calibration.R     # Full factorial calibration
│   ├── benchmarking_helpers.R  # Benchmark targets and comparison functions
│   └── GRAM_calibration_validation_protocol.docx
│
├── analyses/                 # Analysis projects, one subfolder each
│   ├── model_development/    # Model development and paper-1 figures
│   ├── testing_strategies/   # Cognitive testing strategies (active)
│   └── ...
│
├── validation/               # Internal and external validation against cohorts
│
├── archive/                  # Superseded files, kept for the record only
│
└── docs/                     # Extended documentation
    ├── methods_appendix.Rmd  # Parameter definitions, data sources, rationale
    ├── ONBOARDING.md         # Guide for country-specific adaptations
    └── references.bib        # Bibliography for the methods appendix
```

## Key Documentation

| Document | Purpose |
|----------|---------|
| `docs/methods_appendix.Rmd` | Detailed parameter definitions, data sources, and rationale |
| `docs/ONBOARDING.md` | Orientation for new analysts and country-specific adaptations |
| `model/modules/MODULES.md` | Function-level documentation for each simulation module |
| `data/README.md` | Data sources and how each input file was generated |
| `calibration/GRAM_calibration_validation_protocol.docx` | Calibration parameters, targets, goodness-of-fit measure, and validation criteria |
| `analyses/testing_strategies/QUICKSTART_GUIDE.md` | Running and authoring testing-strategy scenarios |
| `archive/README.md` | What each archived file was, and why it is no longer live |

## Model Overview

GRAM simulates individual trajectories through cognitive health states:
- **Healthy** → **MCI** → **Mild/Moderate/Severe Dementia** → **Death**

Each simulated individual has attributes (age, sex, education, APOE4 status, etc.) that influence transition probabilities. The model supports:
- Natural history simulation
- Screening/diagnostic interventions (e.g., Brain Health Assessment)
- Treatment interventions (disease-modifying therapies)
- Cost-effectiveness analysis

## Requirements

`model/setup.R` loads all of these unconditionally, so all are required — sourcing it
will fail if any are missing.

```r
install.packages(c("tidyverse", "scales", "rlang", "flextable", "patchwork"))
```

## Contact

For questions or contributions, contact the project maintainer: sigal.maya@ucsf.edu
