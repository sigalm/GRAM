# ****************************************************************
# ======= GRAM EXTERNAL VALIDATION: FRAMINGHAM HEART STUDY =======
# ****************************************************************

# Run GRAM using the Framingham Heart Study cohort as the starting population.
# Goal is to replicate dementia incidence from: Satizabal 2016 (https://www.nejm.org/doi/full/10.1056/NEJMoa1504327)

# The general model setup is sourced rather than copied, so this script only states what differs
# from it. No microdata is supplied, so the cohort is built from defaults and the parameters below.

source("model/setup.R") # note: this clears the environment, so keep it first
source("model/helpers/source_all.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")

library(tableone)
library(easystats)


######################################## COHORT OVERRIDES ########################################

## Model settings
l.inputs[["n.ind"]] <- 3100                         # number of individuals to simulate
l.inputs[["n.cycle"]] <- 7                          # number of cycles to simulate (2 more than Framingham because last 2 cycles of model flattens out)
l.inputs[["seed_stochastic"]] <- 20250624

## Demographics of the Framingham cohort
l.inputs[["AGE_start_mean"]] <- 72                  
l.inputs[["AGE_start_sd"]] <- 9
l.inputs[["p.SEX_start_male"]] <- 0.44
l.inputs[["p.SEX_start_female"]] <- 0.56
l.inputs[["p.APOE4_start"]] <- c(0.79, 0.21)        # p for non-carrier and carrier, respectively

# Education is reported for the cohort as a whole, not by race, so it is supplied as a marginal
# distribution. Setting p.EDU_start makes it take precedence over the race-specific m.EDU_start.
# Values are years of education standing in for each reported category.
l.inputs[["v.EDU_val"]]   <- c(16, 12, 8)           # college, high school, less than high school
l.inputs[["p.EDU_start"]] <- c(0.63, 0.32, 0.05)    # must add to 1

# race and income not reported, left as default

## Baseline cognitive status: a prevalent cohort, unlike the default all-healthy start
l.inputs[["p.SYN_start"]] <- c(0.82, 0.03, 0.15)    # p for SYN == 0 (normal), 0.5 (TCI) and 1 (impaired)

## Mortality
l.inputs[["p.HCARE_start"]] <- c(0.25, 0.75)
l.inputs[["hr.mort_mci_age"]] <- c(1, 1, 1)
l.inputs[["hr.mort_mod_age"]] <- l.inputs[["hr.mort_sev_age"]] <- c(1, 1, 1)
l.inputs[["m.lifetable"]][nrow(l.inputs[["m.lifetable"]]), ] <- 1   # certain death at the oldest tabulated age

## Progression and incidence, scaled to the Framingham observation period
l.inputs[["r.CDRfast_mean"]] <- 1.6
# 51-element age curve spanning ages 50-100, read at each individual's current age. This cohort
# starts at 70 and runs 7 cycles, so it occupies roughly ages 70-76 of the curve.
l.inputs[["r.CDRslow_mean"]] <- (seq(0, 1, length.out = 51)^1.8660) * (1.8515 * l.inputs[["r.CDRslow_mean"]])
l.inputs[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 1.9869
# (m.hr_mci no longer needs padding with extra rows: f.calc_MCIprob holds the lookup at the
#  oldest tabulated age, which is what the padding did.)


######################################## BENCHMARKS ########################################

# Benchmarks are supplied explicitly rather than taken from the internal-validation defaults,
# so it stays visible which targets this cohort is actually being held to.

benchmarks_framingham <- benchmarks_internal

# Reside time is dropped: 7 cycles is far too short a follow-up to observe a completed spell
# in MCI or dementia, so every duration would be censored.
benchmarks_framingham[["reside_time"]] <- NULL

# NOTE: mortality, prev_by_age and age_of_onset are still the internal US general-population
# targets, not Framingham-specific figures. Replace them as cohort-specific values become
# available. The dementia incidence target from Satizabal 2016 is not compared here -- there
# is no incidence comparison function.


######################################## RUN ########################################

sim_calib <- run_benchmarking(l.inputs, "", benchmarks = benchmarks_framingham)
sim_calib$mort$result_plot
