# ****************************************************************
# ======= GRAM EXTERNAL VALIDATION: FRAMINGHAM HEART STUDY =======
# ****************************************************************

# Run GRAM using the Framingham Heart Study cohort as the starting population.
# Goal is to replicate dementia incidence from: Binder 2019 (https://link.springer.com/article/10.1007/s10654-019-00567-6)
# Note: Used to be Satizabal 2016 (https://www.nejm.org/doi/full/10.1056/NEJMoa1504327), but demoted because of high risk of underestimation.
# Binder uses the same data but analyzes it differently that is likely to be more accurate.

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
l.inputs1[["n.ind"]] <- 10000                        # number of individuals to simulate (reduce sampling noise with higher N)
l.inputs1[["n.cycle"]] <- 8                          # number of cycles to simulate (2 more than intended follow up because last 2 cycles of model is flat and 1st cycle is baseline)

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
l.inputs1[["p.SYN_start"]] <- c(0.91, 0.03, 0.06)    # p for SYN == 0 (normal), 0.5 (TCI) and 1 (impaired)

# Justification for p.SYN_start: 
# 1. Note that parametric GRAM defaults to assigning ALL prevalent impairment to MCI (p.SEV_start is 1,0,0,0)
# 2. Yuan et al 2020 (https://doi.org/10.3233/JAD-20078) Figure 1 gives information on prevalence of MCI in Epoch 4 
# (same sample as Satizabal and Binder). Specifically, if we focus on Window 1:
## a. N(MCI) = 165. Of note, N=346 was assigned MCI or dementia after the window. I include a proportion of them 
## because Epoch 4 in Binder considers all 3 time windows described here. Assuming same distribution of severity, 43% of those 346 
## should be MCI. That gives 148 more people for a total of N(MCI-adjusted) = 313. 
## Denominator is everyone minus dementia. N(dementia) = 106+67+43+(346-148) = 414. N(total) = 5799.
## MCI prevalence is thus 313 / (5799-414) = 5.8%
## I then assign a conservative 3% of TCI. This doesn't really make a huge difference because the majority of dementia cases in a 5-year
## time horizon are due to those already with MCI progressing to dementia (rather than incident cases developing - simply not enough time).
## Remainder start healthy.


######################################## RUN ########################################

sim_ext_val <- f.run(l.inputs1, microdata = NULL, printLevel = 0)


######################################## DEMENTIA INCIDENCE ########################################

library(epiR)
library(survival)

# One row per individual at risk at baseline, with the follow-up each actually contributed. Baseline is
# cycle 1 and the five Framingham follow-up years are cycles 2-6; the 7th cycle is left out.
df_dem <- f.person_time(sim_ext_val, event = "dementia",
                        baseline_cycle = 1, max_followup = 5)

table(df_dem$status)   # 0 = event-free at 5y, 1 = dementia, 2 = died dementia-free


table(status = df_dem$status, followup_yrs = df_dem$time)
sum(df_dem$time)       # person-years, vs. 5 * nrow() if nobody died

# Incidence rate per 1,000 person-years, with an exact CI.
epi.conf(as.matrix(cbind(sum(df_dem$status == 1), sum(df_dem$time))),
         ctype = "inc.rate", method = "exact") * 1000

######################### LIKE-FOR-LIKE COMPARISON WITH THE FHS TARGET #########################

# TARGET: 3.13 per 100 (95% CI 1.58-4.69), five-year cumulative incidence of dementia.
#
# Binder, Balmford & Schumacher (Eur J Epidemiol 2019, doi:10.1007/s10654-019-00567-6)
# reanalysed the FHS data behind Satizabal 2016 in 5118 participants aged 60+, with a
# spline-based penalised-likelihood illness-death multi-state model. Their simulation shows
# Cox models that censor at death or last observation consistently underestimate incidence,
# while the multi-state estimator recovers it; the original Cox figures are therefore biased
# low. Their epoch series is 3.84 / 2.66 / 3.29 / 3.13 per 100 -- essentially flat, and they
# conclude a decline in dementia incidence is not supported. Because the series is flat the
# choice of epoch barely matters; we take Epoch 4 simply for its recency.
#
# Wording note: the paper calls that series "cumulative hazard rates ... per 100 persons", but 
# what they report is actually a cumulative incidence (mentioned in a paranthetical in text).
# In any case, at this magnitude the distinction is negligible because a 5-year
# hazard rate of 3.13 per 100 person over 5 years corresponds to a cumulative incidence of 
# 1 - exp(-0.0313) = 3.08 per 100, a gap of 0.05 against a confidence interval 3 points wide. 
# I treat it here as an incidence (i.e., proportion).

# Death is a competing risk in the target.

# The multi-state target from Binder et al is the right one for GRAM: the model observes true
# state every cycle, with no missed visits and no undiagnosed onset, so it produces true
# incidence by construction and should be held to an estimator that recovers it.
#
# The target is age- and sex-adjusted, and the estimate below is not. Deferred for now: the
# cohort is simulated at the reported FHS mean age and sex mix, so the crude figure is close
# to the adjusted one, but this is an approximation.
#
# Residual mismatch: Because f.module_mortality() runs before the
# health modules, an individual who dies in cycle t can never be recorded as developing
# dementia in cycle t. GRAM structurally excludes the onset-before-death cases the
# multi-state model intended to recover, which means using our data for the incidence calculations
# below likely bias the estimates low. 

# Five-year cumulative incidence, with death as a competing risk:
df_dem$event <- factor(df_dem$status, levels = 0:2,
                         labels = c("event-free", "dementia", "death"))
fit_dem <- survfit(Surv(time, event) ~ 1, data = df_dem)
s5 <- summary(fit_dem, times = 5)

# pstate / std.err column for "dementia" is the five-year cumulative incidence
j    <- which(s5$states == "dementia")
CI5  <- s5$pstate[, j]
CI5_se <- s5$std.err[, j]

c(estimate = CI5, lower = CI5 - 1.96 * CI5_se, upper = CI5 + 1.96 * CI5_se) * 100
                      # per 100 persons, vs. the target 3.13 (1.58-4.69)

# With p(TCI)=0.03 at baseline: 3.23
# If p(TCI)=0:                  2.50
# If p(TCI)=0.06:               3.78

# The two estimators that ignore the competing risk, for contrast: 
1 - exp(-(sum(df_dem$status == 1) / sum(df_dem$time)) * 5)    # constant-hazard
sum(df_dem$status == 1) / nrow(df_dem)                        # simple proportion

