# GRAM Onboarding Guide: Country-Level Adaptations

This document provides everything needed to onboard a new analyst for leading country-level adaptations of the GRAM model.

---

## Part 1: GitHub Workflow for Country Adaptations

### Overview

The workflow uses **feature branches** for each country adaptation. Only the project lead can merge into `main`, ensuring the core model stays stable.

### Roles

| Role | Permissions | Responsibilities |
|------|-------------|------------------|
| **Project Lead (You)** | Admin, merge to `main` | Review PRs, approve merges, maintain core model |
| **Analyst** | Write access | Create branches, push changes, open PRs |

### Setting Up (One-Time, Using GitHub Desktop)

#### For the Project Lead:
1. **Add collaborator**: Go to GitHub.com → Repository → Settings → Collaborators → Add the analyst's GitHub username
2. **Protect main branch** (recommended):
   - Go to Settings → Branches → Add branch protection rule
   - Branch name pattern: `main`
   - Check: "Require a pull request before merging"
   - Check: "Require approvals" (set to 1)
   - This ensures all changes to `main` go through you

#### For the Analyst:
1. **Accept invitation**: Check email for GitHub invitation, click to accept
2. **Clone repository in GitHub Desktop**:
   - File → Clone Repository
   - Select the GRAM repository
   - Choose local path (e.g., `Documents/GitHub/GRAM`)
3. **Open in RStudio**: Open `GRAM.Rproj` from the cloned folder

### Daily Workflow (Using GitHub Desktop)

#### Starting a New Country Adaptation

1. **Sync with main first**:
   - In GitHub Desktop: Click "Fetch origin" to get latest changes
   - If there are changes, click "Pull origin"

2. **Create a new branch**:
   - In GitHub Desktop: Current Branch → New Branch
   - Name it descriptively: `country/brazil`, `country/south-africa`, `country/india`
   - Click "Create Branch"

3. **Work in RStudio**:
   - Make your changes (data files, parameters, analysis scripts)
   - Save files as you work

4. **Commit changes regularly**:
   - In GitHub Desktop: You'll see changed files listed
   - Write a clear commit message (e.g., "Add Brazil mortality data")
   - Click "Commit to country/brazil"

5. **Push to GitHub**:
   - Click "Push origin" to upload your commits

#### Submitting Work for Review

1. **Open a Pull Request (PR)**:
   - In GitHub Desktop: Branch → Create Pull Request
   - This opens GitHub.com in your browser
   - Fill in:
     - **Title**: "Country adaptation: Brazil"
     - **Description**: What you changed and why
   - Click "Create Pull Request"

2. **Wait for review**: The project lead will review and either:
   - Approve and merge
   - Request changes (you'll get an email)

#### If Changes Are Requested

1. Make the requested changes in RStudio
2. Commit and push again (same branch)
3. The PR updates automatically

### Branch Naming Convention

| Type | Pattern | Example |
|------|---------|---------|
| Country adaptation | `country/<country-name>` | `country/brazil` |
| Bug fix | `fix/<description>` | `fix/mortality-calculation` |
| New feature | `feature/<description>` | `feature/add-screening-pathway` |

### What NOT to Do

- **Never commit directly to `main`** — always use a branch
- **Don't modify core model files** (`model/modules/`, `model/simulation.R`) without discussion
- **Don't delete or rename existing data files** — add new ones instead

### Resolving Conflicts

If GitHub Desktop shows a conflict:
1. **Don't panic** — conflicts happen when the same file was edited in two places
2. Click "Open in RStudio" or your text editor
3. Look for conflict markers (`<<<<<<<`, `=======`, `>>>>>>>`)
4. Keep the code you want, delete the markers
5. Save, commit, and push

---

## Part 2: Orientation Session Plan (1.5 Hours)

### Pre-Session Checklist

- [ ] Analyst has GitHub account
- [ ] Analyst has RStudio installed
- [ ] Analyst has GitHub Desktop installed
- [ ] You've added analyst as collaborator on the repository
- [ ] Analyst has cloned the repository

### Session Agenda

#### Introduction (10 min)
- **Goal**: Understand what GRAM does and your role
- Topics:
  - GRAM simulates cognitive impairment progression (Healthy → MCI → Dementia → Death)
  - Country adaptations require changing **data inputs** and **parameters**, not the core model
  - Your workflow: branch → adapt → PR → review → merge

#### Project Structure Tour (15 min)
- **Goal**: Know where everything is

```
GRAM/
├── model/                    # DO NOT MODIFY (core simulation)
│   ├── setup.R               # Parameters you WILL modify (copy first!)
│   ├── simulation.R          # Main engine (don't touch)
│   └── modules/              # Simulation logic (don't touch)
│
├── data/                     # Input data - ADD new country data here
│   ├── mortality/            # Life tables (country-specific)
│   ├── mci_incidence/        # MCI rates (may need country data)
│   └── acs_data/             # Population data (replace for your country)
│
├── analyses/                 # Your analysis scripts go here
│   └── country_<name>/       # Create a folder for each country
│
└── docs/                     # Documentation
```

#### Running the Model (20 min)
- **Goal**: Successfully run a simulation
- Live demo:
  1. Open `GRAM.Rproj` in RStudio
  2. Source the setup files (show the Quick Start from README)
  3. Run a basic simulation
  4. Examine the output structure

```r
# Quick Start
source("model/setup.R")
source("model/helpers/source_all.R")
source("model/simulation.R")

# Run simulation
microdata <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")
results <- f.wrap_run(l.inputs, microdata = microdata)

# Explore results
names(results)
results$fig.progression
```

#### Understanding Parameters (20 min)
- **Goal**: Know which parameters to change for country adaptations
- Walk through `model/setup.R`, which is divided by `##` section headers. Search for the
  header rather than jumping to a line number — the numbers move whenever the file is edited:
  - `## Demographic inputs`: starting age, sex distribution, education, race/ethnicity,
    income, medical burden. Everything between the `Parametric fallback` and
    `End parametric demographic inputs` markers is used **only** when no microdata is
    passed to `f.initialize()`.
  - `## Mortality`: hazard ratios by syndrome and severity, and the life table
  - `## Logistic regression for transition to MCI from Healthy`: MCI incidence and the
    calibrated multipliers `param1`, `param2a`, `param2b`
  - `## Cognitive test scoring and progression`: CDR-SB cutoffs and progression rates
  - `## Cognitive test performance`: test sensitivity and specificity
  - `## Health state utilities`: quality of life weights
  - `## Costs`: country-specific healthcare costs
- Note that `setup.R` holds **uncalibrated** defaults. `calibrate()`, in
  `model/config/calibrate_config.R`, overlays the calibrated values and is what analyses
  actually run with.

#### GitHub Workflow Practice (15 min)
- **Goal**: Comfortable with branching and committing
- Hands-on:
  1. Create a practice branch: `practice/<analyst-name>`
  2. Make a small change (e.g., add a comment to an analysis file)
  3. Commit with a message
  4. Push to GitHub
  5. Open a PR (you can close it without merging)

#### Troubleshooting Common Issues (10 min)
- **Goal**: Know how to debug
- Cover:
  - "Object not found" → Did you source all files in order?
  - "File not found" → Is your working directory set to GRAM root?
  - Model runs but results look wrong → Check parameter values, compare to US baseline
  - GitHub conflict → See conflict resolution steps above

#### Q&A and Next Steps (10 min)
- Answer questions
- Assign first task: Review existing documentation, run the US baseline model
- Schedule follow-up check-in

---

## Part 3: Country Adaptation Checklist

Use this checklist when adapting GRAM for a new country.

### Phase 1: Setup

- [ ] Create branch: `country/<country-name>`
- [ ] Create analysis folder: `analyses/country_<name>/`
- [ ] Copy `model/setup.R` to your analysis folder as starting point for parameters

### Phase 2: Data Collection

#### Required Data (Must Have)

| Data | Source Suggestions | File to Create/Modify |
|------|-------------------|----------------------|
| **Life tables** (mortality by age/sex) | WHO, national statistics office | `data/mortality/lifetable_<country>.RDS` |
| **Population demographics** | Census, DHS surveys | `data/acs_data/<country>_microdata.RDS` or synthetic |

#### Recommended Data (If Available)

| Data | Source Suggestions | Parameter in `setup.R` |
|------|-------------------|------------------------|
| MCI incidence rates | Published studies | `m.hr_mci` |
| Dementia prevalence | ADI reports, local studies | For validation |
| Healthcare costs | WHO-CHOICE, local data | `c.healthy`, `c.mci`, `c.mil`, etc. |
| Healthcare access rates | DHS, local surveys | `p.HCARE_start` |

### Phase 3: Parameter Adaptation

Create a country-specific setup script (e.g., `analyses/country_brazil/setup_brazil.R`):

```r
# Brazil adaptation
# Start with base US parameters
source("model/setup.R")

# Override country-specific parameters
l.inputs[["scenario"]] <- list(
  title       = "Natural progression model - Brazil",
  description = "Adapted for Brazilian population",
  test = NULL
)

# Demographics (example - replace with real data)
l.inputs[["p.SEX_start_male"]] <- 0.48
l.inputs[["p.SEX_start_female"]] <- 0.52

# Education is measured in YEARS, not ordinal levels. p.EDU_start is an optional
# marginal override that must have one probability per entry of v.EDU_val
# (currently 8: 8, 12, 14, 14, 16, 18, 19, 25 years) and sum to 1. Leaving it NULL
# makes the model draw from the race-stratified matrix m.EDU_start instead.
l.inputs[["v.EDU_val"]]   <- c(...)  # country-specific years-of-schooling levels
l.inputs[["p.EDU_start"]] <- c(...)  # same length as v.EDU_val, sums to 1

# Mortality - load country-specific life table
# IMPORTANT: row 1 must be the probability at age 50, with one row per year of age
# thereafter (see the warning below).
l.inputs[["m.lifetable"]] <- as.matrix(
  readRDS("data/mortality/lifetable_brazil.RDS")[, c("m_prob", "f_prob")]
)

# Costs (in local currency or USD)
l.inputs[["c.mci"]] <- 8000  # adjust to local costs
l.inputs[["c.mil"]] <- 15000
# ... etc
```

> **The life table must start at age 50.** `f.update_ALIVE()` looks up mortality with
> `f.age_index()`, which converts an age to a row as `round(AGE) - 50 + 1`. Row 1 is
> therefore assumed to be age 50, row 2 age 51, and so on, one row per single year of age.
> A table that starts at 51, or that uses 5-year age bands, will silently return the wrong
> probability for every individual rather than raising an error. `f.initialize()` also
> derives the oldest attainable age from the table's height
> (`AGE_max <- 50 + nrow(m.lifetable) - 1`), so a short table quietly caps the cohort's age.
> If your source life table starts at birth, subset it to ages 50+ before use; if it is
> banded, expand it to single years first.

### Phase 4: Validation

- [ ] Run baseline (no intervention) scenario
- [ ] Compare dementia prevalence to published estimates for your country
- [ ] Check mortality patterns match life tables
- [ ] Verify age distribution over time looks reasonable
- [ ] Confirm the life table is aligned as described above — plot modelled vs table
      mortality by single year of age and look for a systematic one-year shift

The calibration and validation approach used for the US model, including the targets,
goodness-of-fit measure and acceptance criteria, is described in
`docs/02_calibration_validation_supplement.Rmd`. Country adaptations that
re-calibrate should follow the same structure; `calibration/run_calibration.R` is the
reference implementation.

### Phase 5: Documentation

- [ ] Document all data sources in your analysis folder README
- [ ] Note any assumptions or adaptations made
- [ ] Create a comparison table: US vs. your country parameters

### Phase 6: Submit for Review

- [ ] Push all changes to your branch
- [ ] Open a Pull Request with:
  - Summary of adaptation
  - Data sources used
  - Any issues or questions
  - Validation results

---

## Part 4: Quick Reference

### Key Files

| File | Purpose | Modify? |
|------|---------|---------|
| `model/setup.R` | All model parameters | Copy, then modify copy |
| `model/simulation.R` | Simulation engine | NO |
| `model/modules/*.R` | Module functions | NO |
| `model/helpers/source_all.R` | Loads all code | NO |
| `data/README.md` | Data documentation | Reference only |

### Essential R Commands

```r
# Set working directory (if not using .Rproj)
setwd("/path/to/GRAM")

# Load everything
source("model/setup.R")
source("model/helpers/source_all.R")
source("model/simulation.R")

# Run simulation
results <- f.wrap_run(l.inputs, microdata = microdata)

# View results
results$fig.progression                    # Progression plot
results$aggregated_results_totpop          # Summary statistics

# Check a specific parameter
l.inputs[["m.lifetable"]]
l.inputs[["p.EDU_start"]]
```

### Getting Help

1. **Check documentation first**: `README.md`, `data/README.md`, `model/modules/MODULES.md`
2. **Search the codebase**: Use RStudio's Find in Files (Ctrl/Cmd + Shift + F)
3. **Ask the project lead**: For conceptual questions or if stuck

### Common Parameter Changes for Country Adaptations

| What to Change | Parameter(s) | Notes |
|----------------|--------------|-------|
| Starting age | `AGE_start_mean`, `AGE_start_sd` | |
| Sex distribution | `p.SEX_start_male`, `p.SEX_start_female` | Must sum to 1 |
| Education levels | `v.EDU_val`, `m.EDU_start`, `p.EDU_start` | Education is in **years**. By default drawn from `m.EDU_start` (one column per race/ethnicity level, each summing to 1). Set `p.EDU_start` instead — same length as `v.EDU_val`, summing to 1 — to draw independently of race/ethnicity; it takes precedence when non-`NULL`. |
| Race/ethnicity | `p.RACEETH_start`, `v.RACEETH_val` | Must stay the same length as each other and as the number of columns in `m.EDU_start` |
| Income categories | `p.INCOME_start` | Adjust thresholds in data generation |
| Mortality | `m.lifetable` | Country-specific file; **row 1 must be age 50**, one row per year (see Phase 3) |
| MCI incidence | `m.hr_mci` | Age-specific rates, indexed by age from 50 the same way as the life table |
| Healthcare access | `p.HCARE_start` | |
| All costs | `c.*` parameters | Convert to consistent currency |

---

## Appendix: Example Country Adaptation Structure

```
analyses/
└── country_brazil/
    ├── README.md                    # Document your adaptation
    ├── setup_brazil.R               # Country-specific parameters
    ├── 01_prepare_data.R            # Scripts to prepare input data
    ├── 02_run_scenarios.R           # Run simulations
    ├── 03_analyze_results.R         # Analysis and visualization
    └── results/                     # Output (gitignored)
        ├── figures/
        └── tables/
```

---

*Last updated: July 2026*
*Contact: sigal.maya@ucsf.edu*
