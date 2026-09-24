# ****************************************************************
# ======= 1. GRAM SETUP ======= 
# ****************************************************************

# clear
cat("\014") # clear console
rm(list = ls()) # clear environment
gc() # garbage collection (i.e., clean up memory)
# 
# load necessary packages
library(tidyverse)
library(scales)
library(rlang)
library(flextable)
library(patchwork)
######################################## 1. DEFINE MODEL INPUTS ########################################
l.inputs <- vector(mode = "list", length = 0)

# vector of attribute names (see supplemental 'Attribute names')
l.inputs[["v.attr_names"]] <- c("TIME","ALIVE","AGE","SEX","EDU","RACEETH","INCOME","MEDBUR","APOE4","HCARE",
                                "DX","TX","TCI","SYN","SELECT","BHA","any_BHA_pos","last_BHA_age",
                                "CDR_track","CDR", "MEMLOSS","SEV",
                                "CDRfast_sd1","CDRslow_sd1","CDR_obs","SEV_obs",
                                "PCP","any_PCP_pos","ACCEPT","NP",  # ACCEPT: the slot is here for its a.random column,
                                                                    # the acceptance draw in f.update_BHA; a.out holds NA
                                                                    # (SELECT and BHA already record who declined). Was the
                                                                    # unbuilt PET slot, renamed in place because adding or
                                                                    # removing an attribute re-maps a.random.
                                "TX2","LTC","QALY","COST_test","COST_fu","COST_tx2","COST_care","COST_tx")
l.inputs[["n.attr"]] <- length(l.inputs[["v.attr_names"]])    # number of attributes

# define or describe possible attribute values
l.inputs[["v.ALIVE_val"]]   <- c(0,1)      # 0 = dead, 1 = alive
l.inputs[["v.SEX_val"]]     <- c(1,2)      # 1 = male, 2 = female
l.inputs[["v.EDU_val"]]     <- c(8,12,14,14,16,18,19,25)  # years of education (matches the EDUC coding of the ACS microdata)
l.inputs[["v.RACEETH_val"]] <- c(0,1,2)    # 0 = White, 1 = Black, 2 = Hispanic
l.inputs[["v.INCOME_val"]]  <- c(0,1,2)    # 0 = low (<$9000/y), 1 = medium ($9000-$36000/y), 2 = high (>$36000/y)
l.inputs[["v.MEDBUR_val"]]  <- 0:15        # count of additional health conditions (15 considered, see documentation)
l.inputs[["v.APOE4_val"]]   <- c(0,1)      # 0 = non-carrier, 1 = carrier (hetero- or homozygous)
l.inputs[["v.HCARE_val"]]   <- c(0,1)      # 0 = no regular healthcare provider, 1 = has regular healthcare provider
l.inputs[["v.DX_val"]]      <- c(0,1)      # 0 = no diagnosis, 1 = diagnosed with cognitive impairment (given true impairment)
l.inputs[["v.SYN_val"]]     <- c(0,0.5,1)  # 0 = healthy, 0.5 = transitional cognitive impairment (TCI), 1 = cognitively impaired
l.inputs[["v.SEV_val"]]     <- c(0,1,2,3)  # 0 = MCI, 1 = mild dementia, 2 = moderate dementia, 3 = severe dementia
l.inputs[["v.MEMLOSS_val"]] <- c(0,1)      # 0 = no memory loss, 1 = memory loss
l.inputs[["v.SELECT_val"]]  <- c(0,1)      # 0 = no subjective cognitive concerns, 1 = has subjective cognitive concerns
l.inputs[["v.TX_val"]]      <- c(0,1)      # 0 = Tx off / not provided / stopped, 1 = Tx on / provided / active (DISEASE-MODIFYING)
l.inputs[["v.NP_val"]]      <- c(0,1)      # 0 = negative neuropsych assessment, 1 = positive neuropsych assessment
l.inputs[["v.TX2_val"]]     <- c(0,1)      # 0 = no treatment (non-DMT), 1 = given treatment (non-DMT)
l.inputs[["v.LTC_val"]]  <- c(0,1)      # 0 = not institutionalized / not in long-term care), 1 = institutionalized / in long-term care



######################################## 1.1. USER-DEFINED MODEL SETTINGS ########################################

# model settings

l.inputs[["n.ind"]] <- 10000                               # number of individuals to simulate
l.inputs[["n.cycle"]] <- 50                                # number of cycles to simulate
l.inputs[["seed_stochastic"]] <- 20250624                  # seed for generating random values that drive stochastic parameters
l.inputs[["seed_pa"]] <- 20241022                          # seed for generating random values that drive probabilistic analysis (currently not in use)
l.inputs[["n.psa"]] <- 10                                  # number of PSA iterations (currently not in use)
l.inputs[["r.discount_QALY"]] <- 0.03
l.inputs[["r.discount_COST"]] <- 0.03


######################################## 1.2. EXTERNAL MODEL INPUTS ########################################

## Demographic inputs
# ---- Parametric fallback: only used if microdata not provided to f.initialize() ----
l.inputs[["AGE_start_mean"]] <- 50
l.inputs[["AGE_start_sd"]] <- 0
l.inputs[["p.SEX_start_male"]] <- 0.49
l.inputs[["p.SEX_start_female"]] <- 0.51
l.inputs[["p.RACEETH_start"]] <- c(0.64, 0.14, 0.22)      # p for RACEETH = 0 (white), RACEETH = 1 (Black), and RACEETH = 2 (Hisp)
l.inputs[["p.INCOME_start"]] <- c(0.05, 0.17, 0.78) # p for low, medium, high income, respectively

# Education, in years (v.EDU_val). Supplied either as a race-specific matrix (m.EDU_start,
# one column per level of v.RACEETH_val, each column summing to 1) or, when a cohort has no
# race-specific data, as a single marginal distribution (p.EDU_start). If p.EDU_start is
# non-NULL it takes precedence and education is drawn independently of RACEETH.
# Source: https://www.equityinhighered.org/indicators/u-s-population-trends-and-educational-attainment/educational-attainment-by-race-and-ethnicity/
l.inputs[["m.EDU_start"]] <- matrix(
  c(0.048, 0.095, 0.248,   #  8 years
    0.274, 0.335, 0.327,   # 12
    0.149, 0.181, 0.130,   # 14
    0.111, 0.110, 0.086,   # 14
    0.261, 0.173, 0.145,   # 16
    0.117, 0.081, 0.047,   # 18
    0.017, 0.010, 0.009,   # 19
    0.023, 0.015, 0.008),  # 25
  nrow = 8, byrow = TRUE,
  dimnames = list(NULL, c("white", "black", "hisp")))
l.inputs[["p.EDU_start"]] <- NULL   # optional marginal override; must match length(v.EDU_val)

# Medical burden at baseline. Prevalence of 0-10 additional conditions, and the relative risk
# of having 2+ conditions by education, applied to p.MEDBUR_2plus_start.
# Source: https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-023-15768-8 figure 1
l.inputs[["p.MEDBUR_start_male"]]   <- c(0.155, 0.205, 0.220, 0.175, 0.105, 0.060, 0.045, 0.020, 0.018, 0.005, 0.001)
l.inputs[["p.MEDBUR_start_female"]] <- c(0.175, 0.205, 0.190, 0.175, 0.140, 0.075, 0.040, 0.018, 0.015, 0.002, 0.002)
l.inputs[["p.MEDBUR_2plus_start"]]  <- 0.531
l.inputs[["rr.MEDBUR_2plus_EDU"]]   <- c("college" = 1, "highschool" = 1.32, "lesshighschool" = 1.58)
# ---- End parametric demographic inputs ----

l.inputs[["p.APOE4_start"]] <- c(0.75, 0.25)        # p for non-carrier and carrier, respectively
l.inputs[["p.HCARE_start"]] <- c(0.25, 0.75)        # p for no regular provider and has regular provider, respectively (assumed)
l.inputs[["p.DX_start"]]    <- c(1,0)               # Everyone undiagnosed at start (assumed)

l.inputs[["p.SYN_start"]] <- c(1, 0, 0)             # p for SYN == 0 (normal), SYN == 0.5 (TCI) and SYN == 1 (impaired), respectively
l.inputs[["p.TCI_start"]] <- c(0.5, 0.5)            # for prevalent TCI cases, p of being in the 1st vs 2nd year of the 2-cycle TCI tunnel
l.inputs[["p.MEMLOSS_start"]] <- c(1, 0)            # p for MEMLOSS == 0 (no memloss) and MEMLOSS == 1 (memloss), respectively. Applies to prevalent MCI cases only.
l.inputs[["p.SEV_start"]] <- c(1, 0, 0, 0)          # p for MCI, mild dem, moderate dem, severe dem, respectively. Conditional on SYN == 1.

l.inputs[["p.MEMLOSS_new"]] <- 0.09   # prob of being non-progressive memory loss for new cognitive impairment

l.inputs[["coef_MEDBUR"]] <- 0.1       # 0.2
l.inputs[["amplification_MEDBUR"]] <- 0.025


## Mortality
l.inputs[["hr.mort_mci"]] <- 1.82
l.inputs[["hr.mort_mil"]] <- 2.92
l.inputs[["hr.mort_mod"]] <- 3.85
l.inputs[["hr.mort_sev"]] <- 9.52

l.inputs[["hr.mort_mci_age"]] <- c(1,1,1)
l.inputs[["hr.mort_mil_age"]] <- c(1,1,1)
l.inputs[["hr.mort_mod_age"]] <- c(1,1,1)
l.inputs[["hr.mort_sev_age"]] <- c(1,1,1)


l.inputs[["m.lifetable"]] <- as.matrix(readRDS("data/mortality/non_dementia_mortality_prob_bysex_byage.RDS")[ , c("m_prob_non_dementia", "f_prob_non_dementia")])



## Logistic regression for transition to MCI from Healthy
l.inputs[["m.hr_mci"]] <- array(data = readRDS("data/mci_incidence/mci_incidence_rate_by_age.RDS")[ , 2], dim = c(51,1),
                                dimnames = list(50:100, "r")) / 1000 # divide by 1000 to scale from 1000 person-years to annual rate
# Uncalibrated baselines for the three calibrated parameters. These are the values that make
# the model reproduce its published inputs unchanged: param1 = 1 leaves the Gillis incidence
# hazards as published, and param2a = 1 with param2b = 2 gives a linear age curve whose average
# equals the published mean of 0.6 CDR-SB/year. They are overwritten at the end of this script
# by the calibrated values; they are kept here to document what calibration is departing from.
l.inputs[["param1"]]  <- 1   # multiplier on the age-specific MCI incidence hazards
l.inputs[["param2a"]] <- 1   # curvature of the CDR-SB progression age curve during MCI
l.inputs[["param2b"]] <- 2   # scales the maximum CDR-SB progression rate during MCI

l.inputs[["log_EDU"]] <- log(0.95)
l.inputs[["log_SEX"]] <- log(1)
l.inputs[["log_RACEETHblack"]] <- log(1)
l.inputs[["log_RACEETHhisp"]] <- log(1)
l.inputs[["log_APOE4"]] <- log(1.18)
l.inputs[["log_MEDBUR"]] <- log(1.09)
l.inputs[["log_INCOMEmed"]] <- log(0.80)
l.inputs[["log_INCOMEhi"]] <- log(0.73)

## Cognitive test scoring and progression
l.inputs[["cutoff_CDR"]] <- c("healthy" = 0, 
                              "mci" = 0.5, 
                              "mild" = 4.5, 
                              "moderate" = 9.5,
                              "severe" = 16.5,
                              "max" = 18.0)
# Lower bounds of each severity band. CDR is continuous, so the bands are half-open and a score
# landing exactly on a cutoff belongs to the more severe band: [0,0.5) healthy, [0.5,4.5) MCI,
# [4.5,9.5) mild, [9.5,16.5) moderate, [16.5,18.0] severe. f.update_SEV() gets this by assigning
# bands in increasing order of severity, so the later assignment wins at each boundary.

l.inputs[["r.CDRfast_mean"]] <- 1.6
l.inputs[["r.CDRfast_sd1"]] <- 2.2/sqrt(160)           # individual variation from mean (fast)
l.inputs[["r.CDRslow_mean"]] <- 0.6
l.inputs[["r.CDRslow_sd1"]] <- 1.2/sqrt(358)           # individual variation from mean (slow)
l.inputs[["r.CDR_sd2"]] <- 0                           # observation-level variation in personal trend
# CDR-SB rates above are for people with MCI and/or Alzheimer's disease, from https://pmc.ncbi.nlm.nih.gov/articles/PMC2809036/ table 4

## Health state utilities
# from Table 2 (community-dwelling columns) of https://journals.sagepub.com/doi/full/10.1177/13872877251350381
l.inputs[["u.healthy"]] <- 0.85 
l.inputs[["u.mci"]] <- 0.77
l.inputs[["u.mil"]] <- 0.68
l.inputs[["u.mod"]] <- 0.49
l.inputs[["u.sev"]] <- 0.22

## Costs
l.inputs[["c.healthy"]] <- 0 # direct cost of healthy (annual)
l.inputs[["c.mci"]] <- 13364 # direct cost of MCI (annual; medical + care)
l.inputs[["c.mil"]] <- 26727 # direct cost of mild dementia (annual; medical + care)
l.inputs[["c.mod"]] <- 31644 # direct cost moderate dementia (annual; medical + care)
l.inputs[["c.sev"]] <- 40645 # direct cost of severe dementia (annual; medical + care)


## Scenarios
# Test performance and testing-pathway costs live in model/test_properties.R.
# Treatment parameters live in model/intervention_properties.R.
# Everything that defines a testing STRATEGY -- probs_select, rr.select_prior, the PCP
# follow-up probabilities, ages, intervals and stop rules -- lives in the scenario
# config and must not be given a default here.

l.inputs[["scenario"]] <- list(
  title       = "Natural progression model - US",
  description = "Natural progression of cognitive impairment, no intervention, US",
  test = NULL
)


######################################## 1.3. CALIBRATED PARAMETERS ########################################

# Applied last so they override the uncalibrated baselines set above. The model is therefore
# calibrated by default: no call is needed at the analysis level, which is what the old
# calibrate() function existed to do and what analyses could silently forget to do.
#
# To run uncalibrated, override l.inputs[["param1"]] / [["param2a"]] / [["param2b"]] after
# sourcing this file. A missing file is an error rather than a silent fall-back to the
# uncalibrated baselines, because a run that is quietly uncalibrated is the exact failure
# this arrangement is meant to remove.

if (!file.exists("model/config/calibrated_params.R")) {
  stop("Missing model/config/calibrated_params.R.\n",
       "  It is generated by calibration/write_calibrated_params.R and should be committed.\n",
       "  Restore it from git, or regenerate it from a saved calibration results file.")
}
source("model/config/calibrated_params.R")
l.inputs[names(l.calibrated)] <- l.calibrated

local({
  p <- attr(l.calibrated, "provenance")
  message(sprintf("GRAM calibrated params: param1 = %s, param2a = %s, param2b = %s",
                  l.calibrated$param1, l.calibrated$param2a, l.calibrated$param2b))
  message(sprintf("  run_id %s | commit %s | calibrated %s | GOF %s",
                  p$run_id, p$git_commit, p$calibrated_at, p$gof))
})
