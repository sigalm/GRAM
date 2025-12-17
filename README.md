# GRAM: GBHI Resource Allocation Model

A microsimulation model for cognitive impairment progression and intervention evaluation in aging populations.

## Quick Start

```r
# 1. Open GRAM.Rproj in RStudio (sets working directory automatically)
#    Or manually: setwd("<path-to-GRAM>")

# 2. Source core scripts (in order)
source("model/setup.R")              # Loads l.inputs with parameters
source("model/helpers/source_all.R") # Loads helper functions
source("model/simulation.R")         # Loads simulation functions

# 3. Run a simulation
microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")  # Load microdata
results <- f.wrap_run(l.inputs, microdata = microdata)            # Run simulation
```

## Directory Structure

```
GRAM/
├── model/                    # Core model code
│   ├── setup.R               # Model initialization and parameters (l.inputs)
│   ├── simulation.R          # Main simulation engine
│   ├── modules/              # Core simulation functions
│   │   └── MODULES.md        # Module documentation
│   ├── helpers/              # Utility functions (e.g., formatting, plotting)
│   │   └── source_all.R      # Sources all modules and helpers
│   └── config/               # Scenario configurations
│
├── data/                     # Input data files
│   ├── acs_data/             # Population microdata (ACS-based)
│   ├── mortality/            # Life tables and mortality rates
│   ├── mci_incidence/        # MCI incidence rates by age
│   ├── cogcon/               # Cognitive concern probabilities
│   └── README.md             # Data provenance documentation
│
├── calibration/              # Model calibration and benchmarking files
├── analyses/                 # Analysis projects (by paper/presentation)
│   └── ...                   # Subfolders for each project
│
└── docs/                     # Extended documentation
    └── model_parameters.Rmd  # Parameter definitions, data sources, rationale
```

## Key Documentation

| Document | Purpose |
|----------|---------|
| `model/modules/MODULES.md` | Function-level documentation for each simulation module |
| `docs/model_parameters.Rmd` | Detailed parameter definitions, data sources, and rationale |
| `data/README.md` | Data sources and how each input file was generated |

## Model Overview

GRAM simulates individual trajectories through cognitive health states:
- **Healthy** → **MCI** → **Mild/Moderate/Severe Dementia** → **Death**

Each simulated individual has attributes (age, sex, education, APOE4 status, etc.) that influence transition probabilities. The model supports:
- Natural history simulation
- Screening/diagnostic interventions (e.g., Brain Health Assessment)
- Treatment interventions (disease-modifying therapies)
- Cost-effectiveness analysis

## Requirements

```r
# Core (required)
install.packages(c("tidyverse", "scales"))

# Output formatting (optional, for tables/plots)
install.packages(c("flextable", "patchwork"))
```

## Contact

For questions or contributions, contact the project maintainer: sigal.maya@ucsf.edu
