# GRAM Model: Module Documentation

This document provides an overview of each major module in the `modules/` directory. Each module encapsulates a specific aspect of the simulation. Each module consists of one or more functions that update specific attributes, and one wrapper function that calls all the update functions for that module.

---

## module_init.R: Initialization Module
**Purpose:**
- Initializes the simulation state arrays for random numbers and individual attributes.
- Handles both synthetic and microdata-driven population initialization.

**Functions:**
- `f.initialize(l.inputs, microdata)`: Sets up the 3D arrays for simulation state and random numbers. Assigns starting values to all attributes for each individual, based on input parameters or external microdata.

**Inputs:**
- `l.inputs`: Model input list (attributes, distributions, settings)
- `microdata`: Optional data frame for initializing the population (if not provided, a synthetic population is generated from distributions defined in `l.inputs`)

**Outputs:**
- List containing `a.random` (random number array) and `a.out` (state array)

---

## module_mortality.R: Mortality Module
**Purpose:**
- Updates the alive/dead status of each individual at each cycle, based on mortality risk.

**Attribute Updates:**
- `f.update_ALIVE(alive.lag, v.AGE.lag, v.SYN.lag, v.SEV.lag, m.lifetable, hr.mort_mci, hr.mort_mil, hr.mort_mod, hr.mort_sev, hr.mort_mci_age, hr.mort_mil_age, hr.mort_mod_age, hr.mort_sev_age, random_cycle)`
  - Computes death probability for each individual based on age, cognitive state, and life table. Randomly determines survival or death for each individual based on the computed probability. Returns a vector (0=dead, 1=alive) for the current cycle.

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
  - Keeps income fixed (no change over time). Income is categorical, and we assume there are no changes in one's income group after the age of 50 (starting age of the cohor in base model). 
- `f.update_HCARE(v.HCARE.lag, v.AGE, random_cycle)`
  - Updates healthcare access; uses a random draw to allow for new healthcare access at age 65 for those without it, to reflect the introduction of Medicare.

---

## module_true_health.R: True Health Module
**Purpose:**
- Manages cognitive disease progression and related attributes, including multimorbidity, cognitive status, memory loss, and CDR scores (unobserved).

**Attribute Updates:**
- `f.update_MEDBUR(v.MEDBUR.lag, ...)`
  - Updates medical burden (multimorbidity) for each individual.
- `f.update_APOE4(v.APOE4.lag)`
  - Keeps APOE4 status fixed (genetic risk, no change over time).
- `f.update_SYN(v.SYN.lag, ...)`
  - Updates cognitive syndrome status (e.g., normal, transitional cognitive impairment, cognitively impaired) based on risk factors and previous state.
- `f.update_TCI(v.SYN, v.TCI.lag, n.alive)`
  - Counts how long an individual has been in the transitional cognitive impairment (TCI) tunnel, the `SYN == 0.5` state. Returns 0 for anyone not currently in the tunnel, and otherwise counts up from 1 on the cycle of entry. Tracking elapsed time in a dedicated attribute — rather than inferring it from `SYN` two cycles back — is what lets a cohort initialized with *prevalent* TCI cases progress correctly, since those individuals have no cycle *t-2* to look back to.
- `f.update_MEMLOSS(v.MEMLOSS.lag, ...)`
  - Updates memory loss indicator which flags those with potentially reversible, non-progressive cognitive impairment.
- `f.update_CDR_track(v.SEV.lag, n.alive)`
  - Assigns each individual to a CDR-SB progression group: slow (0) for those with true MCI, fast (1) for those with true dementia.
- `f.update_CDR(v.CDR.lag, ...)`
  - Updates Clinical Dementia Rating (CDR) score (unobserved true score).
- `f.update_SEV(v.SEV.lag, ...)`
  - Updates severity of cognitive impairment for those with cognitive impairment (SYN == 1) as MCI or mild-moderate-severe dementia based on standard CDR-SB categorization.

---

## module_medical_record.R: Medical Record Module
**Purpose:**
- Simulates observed health states, cognitive test results, diagnoses, and clinical assessments for each individual.

**Attribute Updates:**
- `f.update_SELECT(v.SELECT.lag, ...)`
  - Draws whether an individual is selected for testing, from the active scenario's `probs_select` (a probability by true cognitive status, built with `f.select_matrix()`). What this represents depends on the strategy: spontaneous cognitive concern in a reactive program, an EHR risk flag in a selective one, a random opt-in in an inclusive one. Individuals with a prior diagnosis are never selected. Selection is necessary but not sufficient for a test: a selected individual is offered one, and may decline it -- see `f.update_BHA()`.
- `f.update_BHA(v.BHA.lag, ...)`
  - Updates Brain Health Assessment (BHA) test result. Whether an eligible individual is tested in a given cycle is governed by the active scenario's `stop_rule()`, which can retire someone from further testing based on their testing history (e.g. a previous positive result, or reaching a maximum age). See `scenario_TEMPLATE_config.R` for the arguments a `stop_rule` receives. Someone who is due and selected is offered a test and takes it up with the scenario's `probs_accept` (by true cognitive status), or with `p.accept_after_decline` if they declined the last offer they had; without `probs_accept` every offer is taken up. BHA is -9 when not due, -8 when due but not tested, 0/1 for a negative/positive result; for -8, `SELECT` gives the reason (0 not selected, 1 declined). A decline restarts `repeat_interval` just as a test does.
- `f.update_CDR_obs(v.CDR_obs.lag, ...)`
  - Updates observed CDR-SB score (may differ from true CDR-SB).
- `f.update_SEV_obs(v.SEV_obs.lag, ...)`
  - Updates observed severity of cognitive impairment (based on observed CDR-SB score).
- `f.update_DX(v.DX.lag, ...)`
  - Updates clinical/pre-existing diagnosis status (e.g., MCI, dementia, normal). This attribute tracks diagnoses that occur outside of the BHA pathway, such that those with DX == 1 are ineligible for future BHA testing. It reflects empirical estimates of underdiagnosis of MCI and dementia.
- `f.update_PCP(v.PCP.lag, ...)`
  - Updates primary care provider assessment of cognitive status, if applicable.
- `f.update_NP(v.NP.lag, ...)`
  - Updates neuropsychological or specialist assessment of cognitive status, if applicable.

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
- Updated state array (`a.out`) with cognitive and health status attributes

---


For more details, see in-code comments and function documentation within each module.
