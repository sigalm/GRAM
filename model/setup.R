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
                                "DX","TX","TCI","SYN","COGCON","BHA","any_BHA_pos","last_BHA_age",
                                "CDR_track","CDR", "MEMLOSS","SEV",
                                "CDRfast_sd1","CDRslow_sd1","CDR_obs","SEV_obs",
                                "PCP","any_PCP_pos","PET","NP",
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
l.inputs[["v.COGCON_val"]]  <- c(0,1)      # 0 = no subjective cognitive concerns, 1 = has subjective cognitive concerns
l.inputs[["v.TX_val"]]      <- c(0,1)      # 0 = Tx off / not provided / stopped, 1 = Tx on / provided / active (DISEASE-MODIFYING)
l.inputs[["v.PET_val"]]     <- c(0,1)      # 0 = negative PET scan, 1 = positive PET scan
l.inputs[["v.NP_val"]]      <- c(0,1)      # 0 = negative neuropsych assessment, 1 = positive neuropsych assessment
l.inputs[["v.TX2_val"]]     <- c(0,1)      # 0 = no treatment (non-DMT), 1 = given treatment (non-DMT)
l.inputs[["v.LTC_val"]]  <- c(0,1)      # 0 = not institutionalized / not in long-term care), 1 = institutionalized / in long-term care



######################################## 1.1. USER-DEFINED MODEL SETTINGS ########################################

# model settings

l.inputs[["n.ind"]] <- 10000                               # number of individuals to simulate
l.inputs[["n.cycle"]] <- 50                                # number of cycles to simulate
l.inputs[["seed_stochastic"]] <- 20250624                  # seed for generating random values that drive stochastic parameters
l.inputs[["strategy"]] <- NA                               # empty parameter to be filled in as part of the strategies
l.inputs[["strategy_strat1"]] <- "control"
l.inputs[["strategy_strat2"]] <- "intervention_dmt"
l.inputs[["Tx"]] <- 0                                      # empty parameter to be filled in as part of the strategies
l.inputs[["Tx_strat1"]] <- 0
l.inputs[["Tx_strat2"]] <- 1
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
l.inputs[["param1"]]  <- 1   # incidence multiplier (set to calibrated value via calibrate())
l.inputs[["param2a"]] <- 1   # CDR progression curvature exponent (set to calibrated value via calibrate())
l.inputs[["param2b"]] <- 2   # CDR max progression multiplier (set to calibrated value via calibrate()); default = 2 so age-averaged rate equals published mean of 0.6 CDR-SB/year

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
l.inputs[["r.CDR_sd3"]] <- 0                           # rater error (inter-rater reliability, will be less reliable in MCI, better in dem)
# CDR-SB rates above are for people with MCI and/or Alzheimer's disease, from https://pmc.ncbi.nlm.nih.gov/articles/PMC2809036/ table 4

## Cognitive test performance
# Source: Elena Tsoy (both CS and GS at the -1.5z cutoff)
l.inputs[["sens_BHACS"]] <- c(0.44, 0.48, 0.66, 0.96)  # sens[1] for prodromal CI, sens[2] for memory loss (assumed), sens[3] for MCI, sens[4] for dem
l.inputs[["spec_BHACS"]] <- 0.93

l.inputs[["sens_BHAGS"]] <- c(0.56, 0.61, 0.84, 0.98)
l.inputs[["spec_BHAGS"]] <- 0.92

# PCP acts as a downstream filter on BHA-positive patients.
# Combined sens = BHA_sens[sev] * p.PCP_confirm_TP[sev]; combined spec = BHA_spec + (1-BHA_spec) * p.PCP_reject_FP
# Indexed by SEV: [1]=MCI, [2]=mild, [3]=moderate, [4]=severe
l.inputs[["p.PCP_confirm_TP"]] <- c(0.75, 0.88, 0.95, 0.98)
l.inputs[["p.PCP_reject_FP"]]  <- 0.65   # P(PCP correctly dismisses a BHA false positive)


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
l.inputs[["c.bha"]] <- 200   # cost of administering the BHA
l.inputs[["c.bhapos"]] <- 2000   # cost of follow up with patient with positive BHA
l.inputs[["c.Tx"]] <- 5000   # cost of DMT
l.inputs[["c.pet"]] <- 500   # cost of administering a PET scan
l.inputs[["c.np"]] <- 1000  # cost of neuropsych assessment
l.inputs[["c.Tx2"]] <- 500   # cost of non-DMT treatment


## Treatments
l.inputs[["rr.Tx_mci"]] <- 0.70
l.inputs[["Tx_t_max"]] <- 3
l.inputs[["p.Tx"]] <- c(0,1) # Probability of DMT ineligible, vs. eligible
l.inputs[["rr.Px_mci"]] <- 1 # Hypothetical -- risk ratio for developing MCI given a prevention intervention (effectiveness of intervention)

## Scenarios

l.inputs[["m.cogcon_reactive"]] <- readRDS("data/cogcon/m.cogcon_reactive.RDS")
l.inputs[["m.cogcon_selective"]] <- readRDS("data/cogcon/m.cogcon_selective.RDS")
l.inputs[["m.cogcon"]] <- l.inputs[["m.cogcon_reactive"]] %>%
  mutate(h = 1, mci = 1, dem = 1)                   # The default model with not consider cognitive concerns (i.e., everyone has concerns)

l.inputs[["rr.cogcon_prior"]] <- 2    # risk ratio for reporting cognitive concerns if concerns were reported in previous cycle (only acts on t-1)

l.inputs[["scenario"]] <- list(
  title       = "Natural progression model - US",
  description = "Natural progression of cognitive impairment, no intervention, US",
  test = NULL
)
