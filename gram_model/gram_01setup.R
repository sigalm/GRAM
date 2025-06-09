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
library(ggpattern)
library(flextable)
library(knitr)
library(tableone)
library(kableExtra)
library(tibble)
library(patchwork)
######################################## 1. DEFINE MODEL INPUTS ########################################
l.inputs <- vector(mode = "list", length = 0)

# vector of attribute names (see supplemental 'Attribute names')
l.inputs[["v.attr_names"]] <- c("TIME","ALIVE","AGE","SEX","EDU","RACEETH","INCOME","MEDBUR","APOE4","HCARE",
                                "DX","TX","TCI","SYN","COGCON","BHA","CDR_track","CDR", "MEMLOSS", "SEV",
                                "CDRfast_sd1","CDRslow_sd1","CDR_obs", "SEV_obs","FUN",
                                "BEH","INSTIT","QALY","COST_care", "COST_tx")
l.inputs[["n.attr"]] <- length(l.inputs[["v.attr_names"]])    # number of attributes

# define or describe possible attribute values
l.inputs[["v.ALIVE_val"]]   <- c(0,1)      # 0 = dead, 1 = alive
l.inputs[["v.SEX_val"]]     <- c(1,2)      # 1 = male, 2 = female
l.inputs[["v.EDU_val"]]     <- c(1,2,3)    # 1 = college or more, 2 = high school or GED, 3 = less than high school 
l.inputs[["v.RACEETH_val"]] <- c(0,1,2)    # 0 = White, 1 = Black, 2 = Hispanic
l.inputs[["v.INCOME_val"]]  <- c(0,1,2)    # 0 = low (<$9000/y), 1 = medium ($9000-$36000/y), 2 = high (>$36000/y)
l.inputs[["v.MEDBUR_val"]]  <- 0:15        # count of additional health conditions (15 considered, see documentation)
l.inputs[["v.APOE4_val"]]   <- c(0,1)      # 0 = non-carrier, 1 = carrier (hetero- or homozygous)
l.inputs[["v.HCARE_val"]]   <- c(0,1)      # 0 = no regular healthcare provider, 1 = has regular healthcare provider
l.inputs[["v.DX_val"]]      <- c(0,1)      # 0 = no diagnosis, 1 = diagnosed with cognitive impairment (given true impairment)
l.inputs[["v.SYN_val"]]     <- c(0,1)      # 0 = healthy, 1 = cognitively impaired
l.inputs[["v.SEV_val"]]     <- c(0,1,2,3)  # 0 = MCI, 1 = mild dementia, 2 = moderate dementia, 3 = severe dementia
l.inputs[["v.MEMLOSS_val"]] <- c(0,1)      # 0 = no memory loss, 1 = memory loss
l.inputs[["v.COGCON_val"]]  <- c(0,1)      # 0 = no subjective cognitive concerns, 1 = has subjective cognitive concerns
l.inputs[["v.AB_val"]]      <- c(0,1)      # 0 = no amyloid beta (normal), 1 = has amyloid beta (AD pathology)
l.inputs[["v.TX_val"]]      <- c(0,1)      # 0 = Tx off / not provided / stopped, 1 = Tx on / provided / active


######################################## 1.1. USER-DEFINED MODEL SETTINGS ########################################

# model settings
l.inputs[["scenario"]] <- "GRAM natural course of disease" # name of the scenario
l.inputs[["n.ind"]] <- 10000                               # number of individuals to simulate
l.inputs[["n.cycle"]] <- 50                                # number of cycles to simulate
l.inputs[["seed_stochastic"]] <- 20240202                  # seed for generating random values that drive stochastic parameters
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
# model inputs for deterministic analysis

## Demographic inputs
l.inputs[["AGE_start_mean"]] <- 50
l.inputs[["AGE_start_sd"]] <- 0

l.inputs[["p.SEX_start_male"]] <- 0.49
l.inputs[["p.SEX_start_female"]] <- 0.51

l.inputs[["p.EDU_start"]] <- c(0.536, 0.362, 0.102) # p for college, high school, less than high school, respectively. must add to 1.
l.inputs[["p.RACEETH_start"]] <- c(0.64, 0.14, 0.22)      # p for RACEETH = 0 (white), RACEETH = 1 (Black), and RACEETH = 2 (Hisp)
l.inputs[["p.INCOME_start"]] <- c(0.05, 0.17, 0.78) # p for low, medium, high income, respectively
l.inputs[["p.APOE4_start"]] <- c(0.75, 0.25)        # p for non-carrier and carrier, respectively
l.inputs[["p.HCARE_start"]] <- c(0.25, 0.75)        # p for no regular provider and has regular provider, respectively (assumed)
l.inputs[["p.DX_start"]]    <- c(1,0)               # Everyone undiagnosed at start (assumed)

l.inputs[["p.SYN_start"]] <- c(1, 0)                # p for SYN == 0 (normal) and SYN == 1 (impaired), respectively
l.inputs[["p.MEMLOSS_start"]] <- c(1, 0)            # p for MEMLOSS == 0 (no memloss) and MEMLOSS == 1 (memloss), respectively
l.inputs[["p.SEV_start"]] <- c(1, 0, 0, 0)          # p for MCI, mild dem, moderate dem, severe dem, respectively

l.inputs[["m.cogcon_spon"]] <- readRDS("gram_data/cogcon/m.cogcon_spon.RDS")
l.inputs[["m.cogcon_elic"]] <- readRDS("gram_data/cogcon/m.cogcon_elic.RDS")
l.inputs[["m.cogcon"]] <- l.inputs[["m.cogcon_spon"]] %>%
  mutate(h = 1, mci = 1, dem = 1)                   # The default model with not consider cognitive concerns (i.e., everyone has concerns)

l.inputs[["p.MEMLOSS_new"]] <- 0.09   # prob of being non-progressive memory loss for new cognitive impairment

l.inputs[["MEDBUR_start"]] <- readRDS("gram_data/medbur/initial_medbur.RDS")

l.inputs[["coef_MEDBUR"]] <- 0.1       # 0.2
l.inputs[["amplification_MEDBUR"]] <- 0.025


## Mortality
l.inputs[["hr.mort_mci"]] <- 1.82
l.inputs[["hr.mort_mil"]] <- 2.92
l.inputs[["hr.mort_mod"]] <- 3.85
l.inputs[["hr.mort_sev"]] <- 9.52

l.inputs[["hr.mort_mci_age"]] <- c(1,1,1)
l.inputs[["hr.mort_mil_age"]] <- c(1,1,1)
l.inputs[["hr.mort_mod_age"]] <- c(1,0.6,0.3)
l.inputs[["hr.mort_sev_age"]] <- c(1,0.6,0.3)


l.inputs[["m.lifetable"]] <- array(data = as.matrix(
  readRDS("gram_data/mortality/non_dementia_mortality_prob_by_age_v2.RDS")[ , "prob_non_dementia"]), 
  dim = c(51,1), dimnames = list(50:100, "q"))

## Logistic regression for transition to MCI from Healthy
l.inputs[["m.hr_mci"]] <- array(data = readRDS("gram_data/mci_incidence/mci_incidence_rate_by_age.RDS")[ , 2], dim = c(51,1),
                                dimnames = list(50:100, "r")) / 1000 # divide by 1000 to scale from 1000 person-years to annual rate
l.inputs[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 4    # adjust baseline for risk factors

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
                              "max" = 18.0)  # 0 is healthy, 0.5-4.0 is MCI, 4.5-9.0 is mild, 9.5-16.0 mod, 16.5-18.0 severe 

l.inputs[["r.CDRfast_mean"]] <- 1.6
l.inputs[["r.CDRfast_sd1"]] <- 2.2/sqrt(160)           # individual variation from mean (fast)
l.inputs[["r.CDRslow_mean"]] <- 0.6
l.inputs[["r.CDRslow_sd1"]] <- 1.2/sqrt(358)           # individual variation from mean (slow)
l.inputs[["r.CDR_sd2"]] <- 0                           # observation-level variation in personal trend
l.inputs[["r.CDR_sd3"]] <- 0                           # rater error (inter-rater reliability, will be less reliable in MCI, better in dem)
# CDR-SB rates above are for people with MCI and/or Alzheimer's disease, from https://pmc.ncbi.nlm.nih.gov/articles/PMC2809036/ table 4

## Cognitive test performance
# Source: Possin 2018, MCI due to AD vs. control
l.inputs[["sens_BHA"]] <- c(0.50, 0.54, 0.72, 1.00)  # sens[1] for prodromal CI, sens[2] for memory loss (assumed), sens[3] for MCI, sens[4] for dem
l.inputs[["spec_BHA"]] <- 0.85


## Health state utilities
l.inputs[["u.healthy"]] <- 0.85 # PLACEHOLDER - do we want to apply age-related decreases to this?
l.inputs[["u.mci"]] <- 0.73
l.inputs[["u.mil"]] <- 0.69
l.inputs[["u.mod"]] <- 0.53
l.inputs[["u.sev"]] <- 0.38

## Costs
l.inputs[["c.healthy"]] <- 0 # PLACEHOLDER - do we want to apply age-related decreases to this?
l.inputs[["c.mci"]] <- 13364
l.inputs[["c.mil"]] <- 26727
l.inputs[["c.mod"]] <- 31644
l.inputs[["c.sev"]] <- 40645
l.inputs[["c.Tx"]] <- 5000

## Treatments
l.inputs[["rr.Tx_mci"]] <- 0.70
l.inputs[["Tx_t_max"]] <- 3
l.inputs[["p.Tx"]] <- c(0,1) # Probability of DMT ineligible, vs. eligible
l.inputs[["rr.Px_mci"]] <- 1 # Hypothetical -- risk ratio for developing MCI given a prevention intervention (effectiveness of intervention)

