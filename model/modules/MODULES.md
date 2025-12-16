# GRAM Model: Module Documentation

This document provides an overview of each major module in the `gram_helpers/` directory. Each module encapsulates a specific aspect of the simulation, supporting code clarity, modularity, and ease of extension. Each module consists of one or more functions that update specific attributes, and one wrapper function that calls all the update functions for that module.

---

## module_init.R: Initialization Module
**Purpose:**
- Initializes the simulation state arrays for random numbers and individual attributes.
- Handles both synthetic and microdata-driven population initialization.

**Functions:**
- `f.initialize(l.inputs, microdata)`: Sets up the 3D arrays for simulation state and random numbers. Assigns starting values to all attributes for each individual, based on input parameters or microdata.

**Inputs:**
- `l.inputs`: Model input list (attributes, distributions, settings)
- `microdata`: Optional data frame for initializing the population

**Outputs:**
- List containing `a.random` (random number array) and `a.out` (state array)

---

## module_mortality.R: Mortality Module
**Purpose:**
- Updates the alive/dead status of each individual at each cycle, based on mortality risk.

**Attribute Updates:**
- `f.update_ALIVE(alive.lag, v.AGE.lag, v.SYN.lag, v.SEV.lag, m.lifetable, hr.mort_mci, hr.mort_mil, hr.mort_mod, hr.mort_sev, hr.mort_mci_age, hr.mort_mil_age, hr.mort_mod_age, hr.mort_sev_age, random_cycle)`
  - Computes death probability for each individual based on age, syndrome, severity, and life table. Randomly determines survival or death for each individual based on the computed probability. Returns a vector (0=dead, 1=alive) for the current cycle.

---

## module_socdem.R: Sociodemographics Module
**Purpose:**
- Updates sociodemographic attributes (age, sex, education, race/ethnicity, income, healthcare access) over time.

**Attribute Updates:**
- `f.update_AGE(v.AGE.lag)`
  - Increments age by 1 year per cycle.
- `f.update_SEX(v.SEX.lag)`
  - Keeps sex fixed (no change over time).
- `f.update_EDU(v.EDU.lag)`
  - Keeps education fixed (no change over time).
- `f.update_RACEETH(v.RACEETH.lag)`
  - Keeps race/ethnicity fixed (no change over time).
- `f.update_INCOME(v.INCOME.lag)`
  - Keeps income fixed (no change over time). Income is categorical, and we assume no change in one's income group after the age of 50 (starting age of the cohort).
- `f.update_HCARE(v.HCARE.lag, v.AGE, random_cycle)`
  - Updates healthcare access; allows for new healthcare access at age 65 for those without it, using a random draw, to reflect the introduction of Medicare at age 65.

---

## module_true_health.R: True Health Module
**Purpose:**
- Manages cognitive health progression and related attributes, including multimorbidity, cognitive status, memory loss, and CDR scores (unobserved).

**Attribute Updates:**
- `f.update_MEDBUR(v.MEDBUR.lag, ...)`
  - Updates medical burden (multimorbidity) for each individual.
- `f.update_APOE4(v.APOE4.lag)`
  - Keeps APOE4 status fixed (genetic risk, no change over time).
- `f.update_SYN(v.SYN.lag, ...)`
  - Updates syndrome status (e.g., normal, MCI, dementia) based on risk factors and previous state.
- `f.update_MEMLOSS(v.MEMLOSS.lag, ...)`
  - Updates memory loss indicator which flags those with potentially reversible, non-progressive cognitive impairment.
- `f.update_CDR(v.CDR.lag, ...)`
  - Updates Clinical Dementia Rating (CDR) score (unobserved true score).
- `f.update_SEV(v.SEV.lag, ...)`
  - Updates severity of cognitive impairment.

---

## module_medical_record.R: Medical Record Module
**Purpose:**
- Simulates observed health states, cognitive test results, diagnoses, and clinical assessments for each individual.

**Attribute Updates:**
- `f.update_COGCON(v.COGCON.lag, ...)`
  - Updates observed cognitive concern. For scenarios that do not consider cognitive concerns, this attribute is set to 1 for all individuals.
- `f.update_BHA(v.BHA.lag, ...)`
  - Updates Brain Health Assessment (BHA) test result.
- `f.update_CDR_obs(v.CDR_obs.lag, ...)`
  - Updates observed CDR score (may differ from true CDR).
- `f.update_SEV_obs(v.SEV_obs.lag, ...)`
  - Updates observed severity of cognitive impairment.
- `f.update_DX(v.DX.lag, ...)`
  - Updates diagnosis status (e.g., MCI, dementia, normal). This attribute tracks diagnoses that occur outside of the BHA pathway, such that those with DX == 1 are ineligible for future BHA testing. It reflects empirical estimates of underdiagnosis of MCI and dementia.
- `f.update_PCP(v.PCP.lag, ...)`
  - Updates primary care provider assessment of cognitive status.
- `f.update_PET(v.PET.lag, ...)`
  - Updates PET scan result.
- `f.update_NP(v.NP.lag, ...)`
  - Updates neuropsychological or specialist assessment of cognitive status.

---

## module_treatment.R: Treatment Module
**Purpose:**
- Handles treatment and care status updates, including disease-modifying therapy (DMT), non-DMT treatment, and long-term care.

**Attribute Updates:**
- `f.update_TX(v.TX.lag, ...)`
  - Updates disease-modifying therapy (DMT) status.
- `f.update_TX2(v.TX2.lag, ...)`
  - Updates non-DMT treatment status.
- `f.update_LTC(v.LTC.lag, ...)`
  - Updates long-term care status.

**Inputs:**
- `l.inputs`, `a.out`, `t`, `a.random`, `alive`, `n.alive`

**Outputs:**
- Updated state array with cognitive and health status attributes

---


## Additional Helper Modules

- **output_formatters.R**: Functions for formatting simulation outputs (tables, summaries)
- **plotting_helpers.R**: Functions for generating plots and figures
- **test_performance_helpers.R**: Functions for calculating test performance metrics (sensitivity, specificity, etc.)
- **generic_helpers.R, epi_helpers.R, program_utility_helpers.R**: General-purpose utilities for probability, sampling, and code infrastructure
- **run_wrappers.R**: Wrapper functions to run scenarios and aggregate results
- **source_all.R**: Utility to source all helper scripts at once

For more details, see in-code comments and function documentation within each script.
