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

######################################## 1. DEFINE MODEL INPUTS ########################################
l.inputs <- vector(mode = "list", length = 0)

# vector of attribute names (see supplemental 'Attribute names')
l.inputs[["v.attr_names"]] <- c("TIME","ALIVE","AGE","SEX","EDU","RACEETH","INCOME","MEDBUR","APOE4",
                                "TX","SYN","BHA","CDR_track","CDR",
                                "CDRfast_sd1","CDRslow_sd1","CDR_obs",
                                "SEV","SEV_obs","FUN","AGE_MCI",
                                "BEH","INSTIT","QALY","COST_care", "COST_tx")
l.inputs[["n.attr"]] <- length(l.inputs[["v.attr_names"]])    # number of attributes

# define or describe possible attribute values
l.inputs[["v.ALIVE_val"]] <- c(0,1)    # 0 = dead, 1 = alive
l.inputs[["v.SEX_val"]] <- c(1,2)      # 1 = male, 2 = female
l.inputs[["v.EDU_val"]] <- c(1,2,3)    # 1 = college or more, 2 = high school or GED, 3 = less than high school 
l.inputs[["v.RACEETH_val"]] <- c(0,1)  # 0 = not Hispanic or Black, 1 = Hispanic or Black
l.inputs[["v.INCOME_val"]] <- c(0,1,2) # 0 = low (<$9000/y), 1 = medium ($9000-$36000/y), 2 = high (>$36000/y)
l.inputs[["v.MEDBUR_val"]] <- 0:15     # count of additional health conditions (15 considered, see documentation)
l.inputs[["v.APOE4_val"]] <- c(0,1)    # 0 = non-carrier, 1 = carrier (hetero- or homozygous)
l.inputs[["v.SYN_val"]] <- c(0,1)      # 0 = healthy, 1 = cognitively impaired
l.inputs[["v.SEV_val"]] <- c(0,1,2,3)  # 0 = MCI, 1 = mild dementia, 2 = moderate dementia, 3 = severe dementia
l.inputs[["v.AB_val"]] <- c(0,1)       # 0 = no amyloid beta (normal), 1 = has amyloid beta (AD pathology)
l.inputs[["v.TX_val"]] <- c(0,1)       # 0 = Tx off / not provided / stopped, 1 = Tx on / provided / active


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

l.inputs[["p.SEX_start_male"]] <- 0.5
l.inputs[["p.SEX_start_female"]] <- 0.5

l.inputs[["p.EDU_start"]] <- c(0.419, 0.476, 0.105) # p for college, high school, less than high school, respectively. must add to 1.
l.inputs[["p.RACEETH_start"]] <- c(0.65, 0.35)      # p for RACEETH = 0 (not Hisp/Black) and RACEETH = 1 (Hisp/Black)
l.inputs[["p.INCOME_start"]] <- c(0.05, 0.45, 0.50) # p for low, medium, high income, respectively
l.inputs[["p.APOE4_start"]] <- c(0.75, 0.25)        # p for non-carrier and carrier, respectively


## Mortality
l.inputs[["hr.mort_mci"]] <- 1.82
l.inputs[["hr.mort_mil"]] <- 2.92
l.inputs[["hr.mort_mod"]] <- 3.85
l.inputs[["hr.mort_sev"]] <- 9.52

l.inputs[["hr.mort_age"]] <- c(1,0.6,0.3)

l.inputs[["m.lifetable"]] <- array(data = as.matrix(
  readRDS("gram_data/non_dementia_mortality_prob_by_age_v2.RDS")[ , "prob_non_dementia"]), 
  dim = c(51,1), dimnames = list(50:100, "q"))


## Logistic regression for transition to MCI from Healthy
l.inputs[["m.hr_mci"]] <- array(data = readRDS("gram_data/mci_incidence_rate_by_age.RDS")[ , 2], dim = c(51,1),
                                dimnames = list(50:100, "r")) / 1000 # divide by 1000 to scale from 1000 person-years to annual rate

l.inputs[["log_EDU"]] <- log(0.95)
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
l.inputs[["r.CDRfast_sd1"]] <- 2.2/sqrt(160)
l.inputs[["r.CDRslow_mean"]] <- 0.6
l.inputs[["r.CDRslow_sd1"]] <- 1.2/sqrt(358)
l.inputs[["r.CDR_sd2"]] <- 0
l.inputs[["r.CDR_sd3"]] <- 0
# CDR-SB rates above from https://pmc.ncbi.nlm.nih.gov/articles/PMC2809036/ table 4

## Cognitive test performance
# Source: Possin 2018, MCI due to AD vs. control
l.inputs[["sens_BHA"]] <- c(0.72, 0.99)  # MCI vs. control, dementia vs. control
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

