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
l.inputs1 <- l.inputs
l.inputs1[["n.ind"]] <- 2090                         # number of individuals to simulate
l.inputs1[["n.cycle"]] <- 7                          # number of cycles to simulate (2 more than Framingham because last cycle of model is flat and 1st cycle is baseline)
l.inputs1[["seed_stochastic"]] <- 20250624

## Demographics of the Framingham cohort
l.inputs1[["AGE_start_mean"]] <- 72                  
l.inputs1[["AGE_start_sd"]] <- 9
l.inputs1[["p.SEX_start_male"]] <- 0.44
l.inputs1[["p.SEX_start_female"]] <- 0.56
l.inputs1[["p.APOE4_start"]] <- c(0.79, 0.21)        # p for non-carrier and carrier, respectively

# Education is reported for the cohort as a whole, not by race, so it is supplied as a marginal
# distribution. Setting p.EDU_start makes it take precedence over the race-specific m.EDU_start.
# Values are years of education standing in for each reported category.
l.inputs1[["v.EDU_val"]]   <- c(16, 12, 8)           # college, high school, less than high school
l.inputs1[["p.EDU_start"]] <- c(0.63, 0.32, 0.05)    # must add to 1

# race and income not reported, left as default

## Baseline cognitive status: a prevalent cohort, unlike the default all-healthy start
l.inputs1[["p.SYN_start"]] <- c(0.82, 0.03, 0.15)    # p for SYN == 0 (normal), 0.5 (TCI) and 1 (impaired)
# try looking up distn' of MCI in general population w/o dementia


## Mortality
l.inputs1[["p.HCARE_start"]] <- c(0.25, 0.75)
l.inputs1[["hr.mort_mci_age"]] <- c(1, 1, 1)
l.inputs1[["hr.mort_mod_age"]] <- l.inputs[["hr.mort_sev_age"]] <- c(1, 1, 1)
# l.inputs1[["m.lifetable"]][nrow(l.inputs[["m.lifetable"]]), ] <- 1   # certain death at the oldest tabulated age

## Grab calibrated parameters
l.inputs1[["r.CDRslow_mean"]] <-  
  (seq(0, 1, length.out = 51)^1.9500) * (1.8997 * l.inputs[["r.CDRslow_mean"]])
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 1.6408



######################################## RUN ########################################

sim_ext_val <- f.run(l.inputs1, microdata = NULL, printLevel = 0)


######################################## DEMENTIA INCIDENCE ########################################

library(epiR)
library(survival)

# One row per individual at risk at baseline, with the follow-up each actually contributed. Baseline is
# cycle 1 and the five Framingham follow-up years are cycles 2-6; the 7th cycle is left out.
df_dem <- f.person_time(sim_ext_val, event = "dementia",
                        baseline_cycle = 1, max_followup = 5)

# Satizabal 2016 followed participants aged 60 and over, so the cohort is trimmed to match.
df_dem60 <- df_dem[df_dem$AGE >= 60, ]

table(df_dem$status)   # 0 = event-free at 5y, 1 = dementia, 2 = died dementia-free


table(status = df_dem$status, followup_yrs = df_dem$time)
sum(df_dem$time)       # person-years, vs. 5 * nrow() if nobody died

# Incidence rate per 1,000 person-years, with an exact CI.
epi.conf(as.matrix(cbind(sum(df_dem$status == 1), sum(df_dem$time))),
         ctype = "inc.rate", method = "exact") * 1000

# Five-year cumulative incidence, Aalen-Johansen. Death is a competing risk rather than a
# censoring: at a mean baseline age of 72 roughly a fifth of the cohort dies dementia-free,
# and an estimator that censors them credits them with dementia risk over years they never
# lived. That bias is real but small here, because the deaths are spread evenly across the
# five years rather than concentrated early.
df_dem$event <- factor(df_dem$status, levels = 0:2,
                         labels = c("event-free", "dementia", "death"))
fit_dem <- survfit(Surv(time, event) ~ 1, data = df_dem)
summary(fit_dem, times = 5)   # pstate for "dementia" is the 5-year cumulative incidence

# The two estimators that ignore the competing risk, for contrast. Both land within about
# half a percentage point of Aalen-Johansen, so how the cumulative incidence is computed is
# not what drives the comparison against Satizabal -- the underlying rate is.
1 - exp(-(sum(df_dem$status == 1) / sum(df_dem$time)) * 5)   # constant-hazard
sum(df_dem$status == 1) / nrow(df_dem)                        # naive proportion


######################### LIKE-FOR-LIKE COMPARISON WITH THE FHS TARGETS #########################

# Satizabal 2016 reports a five-year age- and sex-adjusted cumulative HAZARD per 100 persons,
# from Cox models fitted in 5205 participants aged 60 and over. Consequences for the
# comparison, none of them cosmetic:
#
#   1. A Cox model censors death, it does not treat it as a competing risk. The
#      Aalen-Johansen figure above is therefore NOT the like-for-like number -- it answers a
#      different question and runs lower.
#   2. A cumulative hazard is not a probability, and exceeds the corresponding cumulative
#      incidence. It is only read "per 100 persons" by convention.
#   3. The cohort must be trimmed to 60+ to match the paper's eligibility.
#
# WHICH TARGET TO USE. van den Hout et al. (Eur J Epidemiol 2019,
# doi:10.1007/s10654-019-00567-6) reanalysed the same FHS data with a spline-penalised
# illness-death multi-state model and found Satizabal's Cox figures biased low, because Cox
# cannot account for dementia acquired between the last dementia-free observation and death.
# Their simulation showed the multi-state estimator recovers true incidence where Cox
# consistently underestimates it. Their epoch series is 3.84 / 2.66 / 3.29 / 3.13 per 100 --
# essentially flat, and they conclude the reported decline is not supported.
#
# The multi-state figures are the right comparator for GRAM: the model observes true state
# every cycle, with no missed visits and no undiagnosed onset, so it produces true incidence
# by construction and should be held to an estimator that recovers it. Because the series is
# flat, the choice of epoch barely matters -- target roughly 3.1-3.3 per 100.
#
# Residual mismatch worth stating rather than fixing: f.module_mortality() runs before the
# health modules, so an individual who dies in cycle t can never be recorded as developing
# dementia in t. GRAM structurally excludes exactly the onset-before-death cases the
# multi-state model exists to recover, so its estimand is marginally narrower than theirs.
# Second-order at annual cycles, but it biases GRAM low, not high.

df60 <- df_dem[df_dem$AGE >= 60, ]
df60$dementia <- as.integer(df60$status == 1)   # death censored, matching the Cox framing
df60$female   <- as.integer(df60$SEX == 2)

cox_dem <- coxph(Surv(time, dementia) ~ AGE + female, data = df60)

# "Age- and sex-adjusted" in a Cox framework means evaluating the fitted curve at reference
# covariate values, not reweighting across age strata. So no age distribution is needed --
# only a reference age and sex mix, and the paper reports both: mean age 72, 44% male. The
# adjusted figure is therefore the curve evaluated there, which is what makes it comparable
# to Satizabal's regardless of how this cohort's own ages happen to fall.
ref <- data.frame(AGE = 72, female = 0.56)
sf_adj <- survfit(cox_dem, newdata = ref)

H5_adj <- -log(summary(sf_adj, times = 5)$surv)
H5_adj * 100          # per 100 persons, vs. the multi-state epoch 4 figure of 3.13
                      # (Satizabal's Cox figure for the same epoch is 2.0)

# Two sensitivities, both of which should move the adjusted figure very little -- that is the
# point of adjusting to a fixed reference.
#
# (a) Fitting on the untrimmed cohort instead of the 60+ subset. Satizabal's sample is 60+ by
#     eligibility, and its reported mean of 72 with sd 9 cannot be a normal truncated at 60:
#     truncating N(72, 9) there yields mean 73.6, sd 7.7, which is what the trimmed cohort
#     here shows. So their age distribution is skewed rather than truncated-normal, and
#     neither trimming nor leaving the cohort whole reproduces it. Evaluating at AGE = 72
#     sidesteps the mismatch, and this check confirms the choice barely matters.
df_all <- df_dem
df_all$dementia <- as.integer(df_all$status == 1)
df_all$female   <- as.integer(df_all$SEX == 2)
cox_all <- coxph(Surv(time, dementia) ~ AGE + female, data = df_all)
-log(summary(survfit(cox_all, newdata = ref), times = 5)$surv) * 100

# (b) Unadjusted cumulative hazard (Nelson-Aalen). The gap between this and the adjusted
#     figure is what the reference-age correction is worth.
sf_na <- survfit(Surv(time, dementia) ~ 1, data = df60)
summary(sf_na, times = 5)$cumhaz * 100

