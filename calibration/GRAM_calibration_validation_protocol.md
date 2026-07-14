---
output:
  word_document: default
  html_document: default
---
# GRAM Calibration and Validation Protocol

## 1. Overview

This document describes the calibration and validation protocol for GRAM. The protocol covers: (1) model calibration via full factorial analysis, (2) internal validation, and (3) external validation.

---

## 2. Calibration

### 2.1 Calibration Parameters

Three parameters are calibrated simultaneously via full factorial analysis:

**Parameter 1 — Incidence multiplier (`param1`)**
A scalar applied to the age-specific MCI baseline incidence hazard vector (derived from Gillis et al. 2019). Its purpose is to remove the average risk factor profile embedded in the Gillis population, since GRAM applies individual-level hazard ratios for education, income, race/ethnicity, APOE4, and medical burden explicitly.

- Range: **0.2 – 1.0**
- Theoretical justification: Values >1 would imply the modeled population is higher-risk than the Gillis meta-analysis population, which is implausible given GRAM models a general population. Values approaching 0 would produce near-zero incidence. Floor of 0.2 reflects a conservative lower bound on plausible downward adjustment.

**Parameter 2a — Progression curvature (`param2a`)**
We assume MCI progression rate increases with age. This parameter is the exponent in the age-dependent MCI progression rate function:
```
r.CDRslow(age) = seq(0, 1, length.out = 51)^param2a × (param2b × r.CDRslow_mean)
```
This controls the shape of how MCI progression rate increases with age.

- Range: **1 – 3**
- Theoretical justification: param2a = 1 produces a linear increase; param2a < 1 produces a concave curve (progression faster at younger ages), which is biologically implausible for a progressive neurodegenerative disease. param2a > 3 implies almost no progression until very old age, which is difficult to justify clinically.

**Parameter 2b — Max progression multiplier (`param2b`)**
A multiplier applied to the published mean MCI progression rate (0.6 CDR-SB points/year, from published literature). It sets the maximum value of the age-dependent progression curve.

- Range: **1 – 2.67**
- Theoretical justification: the floor of 1 anchors the maximum rate at the published mean (lower values would imply the true max is below observed means). The ceiling of 2.67 corresponds to a maximum rate of 1.6 CDR-SB points/year, which is equal to the mean progression rate of dementia.

### 2.2 Calibration Targets

Two population-level targets are used, chosen from sources not used as direct model inputs:

1. **Prevalence by age**: three data points at ages 65, 75, and 88 (roughly mid-points of 10-year age bands) for MCI prevalence from a meta-analysis independent of model inputs (Bai et al); six data points at ages 67, 72, 77, 82, 87, and 95 (roughly mid-points of 5-year age bands) for dementia prevalence from Manly et al (cohort study). Age-band mid-points are used as a standard approximation; this introduces slight underestimation of band averages when prevalence rises steeply with age.

2. **All-cause mortality by 1-year age bands**: age-specific annual all-cause mortality rates derived from actuarial life tables (which gives probabilities) published by the Social Security Administration.

### 2.3 Goodness-of-Fit Measure

The goodness-of-fit (GOF) measure is the **normalized Weighted Sum of Squared Differences (WSSD)**:

For each data point within a target:

$$\text{contribution} = \frac{(\text{model output} - \text{observed})^2}{SE^2}$$

where SE is the standard error of the published estimate. If only 95% CI is available: SE = (upper − lower) / 3.92.

The WSSD is computed separately for each target, then normalized by the number of data points within that target before summing:

$$\text{Total GOF} = \frac{\text{WSSD}_{\text{prev}}}{n_{\text{prev}}} + \frac{\text{WSSD}_{\text{mort}}}{n_{\text{mort}}}$$

Normalization ensures that mortality (many 1-year age bands) and prevalence (3 data points) contribute equally to the total GOF, regardless of the number of data points per target.

The parameter combination yielding the **lowest total GOF** is selected as the calibrated parameter set.

### 2.4 Full Factorial Design

All combinations of parameter values across their ranges are evaluated. Initial calibration uses **6 steps per parameter** (216 total combinations), parallelized across available CPU cores. If the optimum falls at a grid boundary, the range is expanded and the analysis is repeated. Grid resolution may be increased in a second pass if warranted.

**Computing**: runs parallelized locally using the `parallel` package in R (8 cores available). Each run uses n = 100,000 individuals; if runtime is prohibitive, n = 50,000 is acceptable given stable outcome estimates at that sample size (~4,500 MCI cases at age 70).

---

## 3. Internal Validation

Internal validation assesses whether the calibrated model produces biologically plausible disease dynamics, distinct from the population-level targets used for calibration.

### 3.1 Face Validity Checks

The following outputs depend on calibrated parameters but are informative as face validity checks:

- **Duration in MCI stage**: depends on param2a and param2b
- **Average age of MCI onset**: depends on param1

These are not independent tests of model fit but confirm that calibrated parameter values produce clinically reasonable disease trajectories.

### 3.2 Independent Validation Target

**Duration in each dementia stage** (mild, moderate, severe): the dementia progression rate is fixed from a published paper (mean and SD applied at the individual level) and is not a calibrated parameter. Duration in each dementia stage is therefore independent of the calibrated parameters and constitutes an internal validation test.

### 3.3 Acceptance Criteria

Face validity checks and internal validation targets are evaluated qualitatively, based on expert input, given limited flexibility of model after calibration.

---

## 4. External Validation

External validation tests the model's predictive validity in an independent population not used for calibration or internal validation.

### 4.1 Approach

1. Define a new starting cohort based on the baseline characteristics of a clinical trial control arm
2. Run GRAM with the calibrated parameter set
3. Compare model-predicted endpoints to the observed clinical trial control arm endpoints at the same follow-up duration

### 4.2 Acceptance Criteria

Acceptance thresholds for external validation endpoints are **pre-specified before model runs**. Primary criterion is that the model output falls within the 95% CI of the trial's observed endpoint.
If CI is not available, we will accept absolute difference <10% of the observed value.

---

*Protocol developed June 2026*
