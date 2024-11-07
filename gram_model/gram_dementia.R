# # clear
# cat("\014") # clear console
# rm(list = ls()) # clear environment
# gc() # garbage collection (i.e., clean up memory)
# 
# # technical
# library(tidyverse)
# save_output <- TRUE # save output to csv files on disk yes/no

######################################## 1. DEFINE MODEL INPUTS ########################################

# vector of attribute names (see supplemental 'Attribute names')
v.attr_names <- c("TIME","ALIVE","AGE","SEX","EDU","RACEETH","INCOME","MEDBUR","APOE4","TX","SYN","BHA","CDR_track","CDR",
                  "CDRfast_sd1","CDRslow_sd1","CDR_obs",
                  "SEV","SEV_obs","FUN","AGE_MCI",
                  "BEH","INSTIT","QALY","COST_care", "COST_tx")
n.attr <- length(v.attr_names) # number        of attributes

# define or describe possible attribute values
v.ALIVE_val <- c(0,1)    # 0 = dead, 1 = alive
v.SEX_val <- c(1,2)      # 1 = male, 2 = female
v.EDU_val <- c(1,2,3)    # 1 = college or more, 2 = high school or GED, 3 = less than high school 
v.RACEETH_val <- c(0,1) # 0 = not Hispanic or Black, 1 = Hispanic or Black
v.INCOME_val <- c(0,1,2)  # 0 = low (<$9000/y), 1 = medium ($9000-$36000/y), 2 = high (>$36000/y)
v.MEDBUR_val <- 0:15      # count of additional health conditions (15 considered, see documentation)
v.APOE4_val <- c(0,1)    # 0 = non-carrier, 1 = carrier (hetero- or homozygous)
v.SYN_val <- c(0,1)     # 0 = healthy, 1 = cognitively impaired
v.SEV_val <- c(0,1,2,3)    # 0 = MCI, 1 = mild dementia, 2 = moderate dementia, 3 = severe dementia
v.AB_val <- c(0,1)         # 0 = no amyloid beta (normal), 1 = has amyloid beta (AD pathology)
v.TX_val <- c(0,1)       # 0 = Tx on / provided / active, 1 = Tx off / not provided / stopped


######################################## 1.1. USER-DEFINED MODEL SETTINGS ########################################

# start with empty input list
l.inputs <- vector(mode = "list", length = 0)

# model settings
l.inputs[["scenario"]] <- "GRAM natural course of disease" # name of the scenario
l.inputs[["n.ind"]] <- 3000 # number of individuals to simulate
l.inputs[["n.cycle"]] <- 50 # number of cycles to simulate
l.inputs[["seed_stochastic"]] <- 20240202 # seed for generating random values that drive stochastic parameters (e.g., determines baseline sex of a specific individual, or whether the event of death occurs for a certain individual at a certain time)
l.inputs[["strategy"]] <- NA # empty parameter to be filled in as part of the strategies
l.inputs[["strategy_strat1"]] <- "control"
l.inputs[["strategy_strat2"]] <- "intervention_diseasemodifying"
l.inputs[["Tx"]] <- 0 # empty parameter to be filled in as part of the strategies
l.inputs[["Tx_strat1"]] <- 0
l.inputs[["Tx_strat2"]] <- 1
l.inputs[["seed_pa"]] <- 20241022 # seed for generating random values that drive probabilistic analysis
l.inputs[["n.psa"]] <- 10 # number of PSA iterations
l.inputs[["r.discount_QALY"]] <- 0.03
l.inputs[["r.discount_COST"]] <- 0.03


######################################## 1.2. EXTERNAL MODEL INPUTS ########################################
# model inputs for deterministic analysis

## Demographic inputs
l.inputs[["AGE_start_mean"]] <- 50
l.inputs[["AGE_start_sd"]] <- 0

l.inputs[["p.SEX_start_male"]] <- 0.5
l.inputs[["p.SEX_start_female"]] <- 0.5

l.inputs[["p.EDU_start"]] <- c(0.4, 0.5, 0.1) # p for college, high school, less than high school, respectively. must add to 1.
l.inputs[["p.RACEETH_start"]] <- c(0.65, 0.35) # p for RACEETH = 0 (not Hisp/Black) and RACEETH = 1 (Hisp/Black)
l.inputs[["p.INCOME_start"]] <- c(0.35, 0.55, 0.1) # p for low, medium, high income, respectively
l.inputs[["p.APOE4_start"]] <- c(0.7, 0.3) # p for non-carrier and carrier, respectively


## Mortality
l.inputs[["hr.mort_mci"]] <- 1.82
l.inputs[["hr.mort_mil"]] <- 2.92
l.inputs[["hr.mort_mod"]] <- 3.85
l.inputs[["hr.mort_sev"]] <- 9.52

l.inputs[["m.lifetable"]] <- array(data = read.csv("../gram_data/table_prob_die_next_year.csv")[ , 2], dim = c(51,1),
                                   dimnames = list(50:100, "q"))

## Logistic regression for transition to MCI from Healthy

l.inputs[["m.hr_mci"]] <- array(data = readRDS("../gram_data/mci_incidence_rate_by_age.RDS")[ , 2], dim = c(51,1),
                                dimnames = list(50:100, "r")) / 1000 # Divide by 1000 to scale from 1000 person-years to annual rate

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

l.inputs[["r.CDRfast_mean"]] <- 1.43
l.inputs[["r.CDRfast_sd1"]] <- 0.3
l.inputs[["r.CDRslow_mean"]] <- 1.1
l.inputs[["r.CDRslow_sd1"]] <- 0.4
l.inputs[["r.CDR_sd2"]] <- 1
l.inputs[["r.CDR_sd3"]] <- 1


## Cognitive test performance
# Source: Possin 2018, MCI due to AD vs. control
l.inputs[["sens_BHA"]] <- 0.72
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
l.inputs[["rr.Tx_mci"]] <- 0.8
l.inputs[["Tx_t_max"]] <- 5






######################################## 1.3. ESTIMATED ########################################
# n/a



######################################## 1.4. CALIBRATED ########################################
# n/a



######################################## 1.5. RANDOM VALUES ########################################

# random values
a.random <- array(data = NA, dim = c(l.inputs[["n.cycle"]], n.attr, l.inputs[["n.ind"]]), dimnames = list(NULL, v.attr_names, NULL)) # generate 'empty' array (i.e., containing NA) with rows (1st dimension) = cycles, columns (2nd dimension) = attributes, 3rd dimension = individuals
set.seed(l.inputs[["seed_stochastic"]]) # set the seed using the earlier defined seed
a.random[,,] <- runif(n = length(a.random)) # put random value from uniform distribution in each element of the array (i.e., replace NA with random value)


######################################## 2. DECISION MODEL IMPLEMENTATION ########################################

######################################## 2.1. FUNCTIONS: GENERIC ########################################
# This section defines all generic functions being used through the model.


######################################## convert probability to different time

# formula from https://doi.org/10.1007/s40273-020-00937-z

# p = probability
# t_new = targeted time period
# t_old = current time period
# RR = risk ratio (e.g. relative risk or hazard ratio)

f.adjustprobability <- function(p, t_new, t_old, RR) {
  rate1 <- -log(1-p)                  # convert probability to rate
  rate2 <- rate1 * RR * (t_new/t_old) # adjust rate (time and risk ratio)
  out <- 1-exp(-rate2)                # convert rate to probability
  return(out)
}

######################################## sample from categorical variable using pre-defined uniform random values

# p_rand = random values from uniform distribution
# p_cat = probabilities of categorical variable
# values = values corresponding to the probabilities in p_cat (must be in same order)

f.qcat <- function(p_rand, p_cat, values = NULL) {
  
  # see supplemental 'Function f.qcat' for details
  
  if (length(p_rand)==0 | length(p_cat)==0 | any(is.na(p_rand)) | any(is.na(p_cat))) { return(NULL) # if input is NULL
  } else if (sum(p_cat)==0) { return(0) # if input is zero
  } else if (is.vector(p_cat)) { # if p_cat is of type vector
    # checks
    if( (sum(p_cat)<0.99999999 | sum(p_cat)>1.00000001) & sum(p_cat)!=0 & !anyNA(p_cat) ) { warning("vector: p_cat do not sum up to 1"); print(p_cat) }
    if(anyNA(p_cat)) warning("vector: p_cat contains NA")
    if(!is.null(values) & length(p_cat)!=length(values)) stop("vector: length of p_cat does not match length of values")
    if(length(p_cat)<=1) stop("vector: 2 or more p_cat values are required")
    if(is.null(values)) values <- 1:length(p_cat)
    
    # function
    out <- cut(x = p_rand, breaks = c(0,cumsum(p_cat)), labels = FALSE, include.lowest = TRUE) # categorize random values based on user-defined breaks
    if(!is.null(values)) out <- values[out] # apply user-defined labels
    return(out)
    
  } else if (is.matrix(p_cat)) { # if p_cat is of type matrix
    # checks
    if(!is.matrix(p_cat)) stop("p_cat is not defined as matrix")
    if(any(is.na(colSums(p_cat))) ) warning("some colsums contain NA")
    if(!all(colSums(p_cat)-1<10E-10) ) warning("matrix: p_cat do not sum up to 1, might be due to NA")
    if(!is.null(values) & nrow(p_cat)!=length(values) ) stop("number of rows of p_cat matrix does not match length of values")
    if(is.null(values) ) values <- 1:nrow(p_cat)
    if(nrow(p_cat)<=1 ) stop("2 or more p_cat values are required (i.e. matrix with 2 or more columns")
    if(length(p_rand)!=ncol(p_cat) ) stop("length of p_rand differs from number of columns of p_cat")
    
    # function
    out <- rep(NA, length(p_rand))
    p_cat.cumsum <- apply(X=p_cat, MARGIN=2, FUN=cumsum)
    n <- nrow(p_cat.cumsum)-1
    out[p_rand>=0 & p_rand<=p_cat.cumsum[1,]] <- values[1]
    for (i in 1:n) {
      out[p_rand>p_cat.cumsum[i,] & p_rand<=p_cat.cumsum[i+1,]] <- values[i+1]
    }
    return(out)
    
  } else stop("p_cat is not vector or matrix")
  
}

######################################## weibull

# mu = location parameter
# gamma = shape parameter

f.predict_weibull <- function(time, mu, gamma) {
  exp(-exp(mu) * time^gamma)
}


######################################## discounting

f.discount <- function(x, l.inputs, discount_rate, n.cycle) {
  # x can be scalar, vector or matrix
  as.matrix(x) / (1 + discount_rate)^(0:(n.cycle - 1))
}


######################################## 2.2. FUNCTIONS: ATTRIBUTE UPDATE ########################################
# This section defines the function for each attribute to be updated.


######################################## TIME

f.update_TIME <- function(v.TIME.lag) {
  time <- v.TIME.lag + 1 # time <- 'time at previous cycle' + 'fixed cycle length of 1 year'
  return(time)
}

########## !!!!!!!!!!!!!!! The model is run by updating each attribute in a loop over the cycles (over time). 
# At each cycle attributes are updated using the attribute status at the previous cycle or the status at the current cycle. 
# Except the first cycle, which is manually put in (i.e., starting values). 
# For transparency, no information from other than the previous or current is used. Information from more than 1 cycle ago 
#     could be used by tracking the history of an attribute in a separate attribute. For each attribute a function is written to update it. 
#     Then, in a loop all functions are called to update their status cycle by cycle. 


######################################## ALIVE

# https://adv-r.hadley.nz/subsetting.html vectorized subset from array (i.e. create matrix with array coordinates in each row to be looked up) 

f.update_ALIVE <- function(
    v.AGE.lag, v.SYN.lag, v.SEV.lag, 
    m.lifetable, 
    hr.mort_mci, hr.mort_mil, hr.mort_mod, hr.mort_sev, 
    random_cycle) 
{
  
  # generate life table coordinates for looking up age-specific mortality 
  # (see https://adv-r.hadley.nz/subsetting.html paragraph 4.2.3 subsetting > selecting multiple elements > subsetting)
  lifetable_lookup_coordinates <- matrix(data = round(v.AGE.lag,0), ncol = 1) - 50 + 1
  
  # determine relative mortality risk related to syndrome and severity
  v.risk_ratio <- rep(NA, length(v.AGE.lag))
  v.risk_ratio[v.SYN.lag==0] <- 1
  v.risk_ratio[v.SYN.lag==1] <- hr.mort_mci
  v.risk_ratio[v.SYN.lag==2 & v.SEV.lag==1] <- hr.mort_mil
  v.risk_ratio[v.SYN.lag==2 & v.SEV.lag==2] <- hr.mort_mod
  v.risk_ratio[v.SYN.lag==2 & v.SEV.lag==3] <- hr.mort_sev
  
  # calculate death probability using life table, relative risk and adjustment for cycle length
  prob_death <- f.adjustprobability(
    p = m.lifetable[lifetable_lookup_coordinates], 
    t_new = 1, 
    t_old = 1,
    RR = v.risk_ratio)
  
  # compare probability to random value
  alive <- as.numeric(!prob_death > random_cycle)
  
  # check for errors
  if (is.na(sum(alive))) {
    print(table(alive, useNA = "always"))
    print(table(alive, v.SYN.lag, useNA = "always"))
    print(table(alive, v.SEV.lag, useNA = "always"))
    stop("NAs produced in 'ALIVE'")
  }
  # return
  return(alive)
}


######################################## AGE

f.update_AGE <- function(v.AGE.lag) {
  age <- v.AGE.lag + 1 # 1 represents fixed cycle time of 1 year
  return(age)
}

######################################## SEX

f.update_SEX <- function(v.SEX.lag) {
  sex <- v.SEX.lag
  return(sex)
}

######################################## RACEETH

f.update_RACEETH <- function(v.RACEETH.lag) {
  raceeth <- v.RACEETH.lag
  return(raceeth)
}

######################################## INCOME

f.update_INCOME <- function(v.INCOME.lag) {
  income <- v.INCOME.lag
  return(income)
}

######################################## EDU

f.update_EDU <- function(v.EDU.lag) {
  edu <- v.EDU.lag
  return(edu)
}

######################################## MEDBUR

f.update_MEDBUR <- function(v.MEDBUR.lag) {
  medbur <- v.MEDBUR.lag
  return(medbur)
}

######################################## APOE4

f.update_APOE4 <- function(v.APOE4.lag) {
  apoe4 <- v.APOE4.lag
  return(apoe4)
}

######################################## MCI

f.update_SYN <- function(l.inputs, v.AGE.lag, v.EDU.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag, v.SYN.lag, random_cycle, n.alive) {
  
  # start with empty vector
  symptoms <- rep(NA, n.alive)

  # look up hazard rate given age
  hazards.age <- l.inputs[["m.hr_mci"]][matrix(data = round(v.AGE.lag,0), ncol = 1) - 50 + 1]
  
  # calculate hazard at current cycle given coefficients, convert to probability
  hazard <- hazards.age * exp(
    l.inputs[["log_EDU"]] * v.EDU.lag + 
    l.inputs[["log_APOE4"]] * v.APOE4.lag + 
    l.inputs[["log_MEDBUR"]] * v.MEDBUR.lag + 
    l.inputs[["log_INCOMEmed"]] * (v.INCOME.lag == 1) + 
    l.inputs[["log_INCOMEhi"]] * (v.INCOME.lag == 2))
  
  prob_none_to_mci <- 1 - exp(-hazard)
  
  # apply dependent on previous state
  symptoms[v.SYN.lag == 0] <- as.numeric(prob_none_to_mci[v.SYN.lag==0] > random_cycle[v.SYN.lag==0])
  symptoms[v.SYN.lag == 1] <- 1
  
  return(symptoms)
}

######################################## BHA

f.update_BHA <- function(v.SYN, v.SYN.lag, sens_BHA, spec_BHA, random_cycle, n.alive) {
    
  bha <- rep(NA, n.alive)

  healthy <- v.SYN == 0
  bha[healthy] <- as.numeric((1-spec_BHA) > random_cycle[healthy])  # If healthy, specificity gives the probability of positive BHA
  
  mci_new <- v.SYN == 1 & v.SYN.lag == 0
  bha[mci_new] <- as.numeric(sens_BHA > random_cycle[mci_new])  # If new MCI, sensitivity gives the probability of positive BHA
  
  mci_still <- v.SYN == 1 & v.SYN.lag == 1
  bha[mci_still] <- 1                                           # If already impaired, positive BHA by default (no re-test)
  #TODO: discuss if this is the right approach.
  return(bha)
}

######################################## CDR-SB fast/slow track

f.update_CDR_track <- function(v.MEDBUR, random_cycle, n.alive) {
  
  track <- rep(NA, n.alive)
  
  x <- 0.5 * v.MEDBUR
  prob_fast <- exp(x) / (1 + exp(x))
  
  track <- prob_fast > random_cycle
  
  return(track)
}

######################################## CDR-SB - true disease status

f.update_CDR <- function(v.SYN, v.SYN.lag, cutoff_CDR, v.CDR.lag, 
                              r.CDRfast_mean, r.CDRslow_mean, v.CDR_track,
                              v.CDRfast_sd1, v.CDRslow_sd1, r.CDR_sd2, random_cycle, n.alive) {
  
  # CDR-SB cut-off values from: O'Bryant et al 2012 (PMC3409562) Table 2
  
  cdr <- rep(NA, n.alive)
  
  # Assign zero to those who are healthy
  healthy <- v.SYN == 0
  cdr[healthy] <- qunif(p = random_cycle[healthy], min = cutoff_CDR["healthy"], max = cutoff_CDR["mci"])
  
  # Assign initial CDR-SB score for new cases. Assume all enter in mild cognitive impairment, skewed right
  mci_new <- v.SYN == 1 & v.SYN.lag == 0
  cdr[mci_new] <- qunif(p = random_cycle[mci_new], min = cutoff_CDR["mci"], max = cutoff_CDR["mild"])  # TODO: Skew right.
  
  # Progress CDR-SB score for those already with impairment
  mci_still <- v.SYN == 1 & v.SYN.lag == 1 
  cdr[mci_still] <- v.CDR.lag[mci_still] + 
    (r.CDRfast_mean * v.CDR_track[mci_still] + r.CDRslow_mean * (1-v.CDR_track[mci_still])) +                               # Add mean increase in CDR-SB score
    (v.CDRfast_sd1[mci_still] * v.CDR_track[mci_still] + v.CDRslow_sd1[mci_still] * (1-v.CDR_track[mci_still])) +  # Add individual-level deviation from mean
    qnorm(p = random_cycle[mci_still], mean = 0, sd = r.CDR_sd2)                                    # Add within-individual deviation
  
  return(cdr)
}


######################################## CDR-SB - observed disease status

f.update_CDR_obs <- function(v.BHA, v.CDR, r.CDR_sd3, random_cycle, n.alive) {
  
  cdr_obs <- rep(NA, n.alive)
  
  cdr_obs[v.BHA == 1] <- v.CDR[v.BHA == 1] + 
    qnorm(p = random_cycle[v.BHA == 1], mean = 0, sd = r.CDR_sd3)        # Add measurement error (rater reliability)
  
  return(cdr_obs)
}


######################################## SEV

f.update_SEV <- function(v.SYN, v.CDR, cutoff_CDR, n.alive) {
  
  # Assign dementia severity based on cognitive score
  sev <- rep(NA, n.alive)
  
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["mci"], cutoff_CDR["mild"])] <- 0
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["mild"], cutoff_CDR["moderate"])] <- 1
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["moderate"], cutoff_CDR["severe"])]<- 2
  sev[v.SYN == 1 & v.CDR >= cutoff_CDR["severe"]] <- 3
  
  return(sev)
}

######################################## SEV_obs

f.update_SEV_obs <- function(v.CDR_obs, cutoff_CDR, n.alive) {
  
  # Assign dementia severity based on cognitive score
  sev_obs <- rep(NA, n.alive)
  
  sev_obs[between(v.CDR_obs, cutoff_CDR["mci"], cutoff_CDR["mild"])] <- 0
  sev_obs[between(v.CDR_obs, cutoff_CDR["mild"], cutoff_CDR["moderate"])] <- 1
  sev_obs[between(v.CDR_obs, cutoff_CDR["moderate"], cutoff_CDR["severe"])]<- 2
  sev_obs[v.CDR_obs >= cutoff_CDR["severe"]] <- 3
  
  return(sev_obs)
}

######################################## FAQ

f.update_FAQ <- function(l.inputs, v.SYN.lag, v.bSEV, v.APOE4) {
  
  # Determine intercept and slope based on health state
  if (v.SYN.lag == 0) {
    intercept <- predictFAQ_int_NL + predictFAQ_bSEV_intercept * v.bSEV + predictFAQ_BSV_intercept_NL
    slope <- predictFAQ_slope_NL + predictFAQ_APOE4 * v.APOE4 + predictFAQ_BSV_slope_NL
  } else if (v.SYN.lag == 1) {
    intercept <- predictFAQ_int_AD + predictFAQ_bSEV_intercept * v.bSEV + predictFAQ_BSV_intercept_MCI
    slope <- predictFAQ_slope_AD + predictFAQ_APOE4 * v.APOE4 + predictFAQ_BSV_slope_MCI
  } else if (v.SYN.lag == 2) {
    intercept <- predictFAQ_int_AD + predictFAQ_bSEV_intercept * v.bSEV + predictFAQ_BSV_intercept_AD
    slope <- predictFAQ_slope_AD + predictFAQ_APOE4 * v.APOE4 + predictFAQ_BSV_slope_AD
  } else {
    stop("Invalid syndrome value (0 = normal, 1 = MCI, 2 = AD).")
  }
  
  # Calculate transformed FAQ score for next cycle, then back-transform
  faq_transformed <- intercept + slope * year + rnorm(1, mean = 0, sd = sqrt(predictFAQ_epsilon))
  faq <- faq_transformed^(1/predictFAQ_kTransform)
  
  # return
  return(faq)
  
  # Can also return categorical for those with AD, if we know cut off points.
}

######################################## TX

f.update_TX <- function(v.TX.lag, v.SYN, v.time, n.alive, random_cycle, Tx_t_max, p.Tx_discontinuation, Tx) {
  
  # start with empty vector
  tx <- rep(NA, n.alive)
  
  # previous Tx state * previous syndrome==MCI * time shorter than maximum treatment duration
  if(Tx==0) tx <- 0
  if(Tx==1) tx <- v.TX.lag * as.numeric(v.SYN==1) * as.numeric(v.time<Tx_t_max)
  
  # return
  return(tx)
  
}


######################################## 2.3. FUNCTIONS: RUN MODEL ########################################

# run model
f.run <- function(a.random, l.inputs) {
  
  # progress bar
  ptm <- proc.time()
  stime <- Sys.time()
  pb = txtProgressBar(min = 0, max = l.inputs[["n.cycle"]], initial = 0)
  
  # empty output array
  a.out <- array(data = NA, dim = c(l.inputs[["n.cycle"]], n.attr, l.inputs[["n.ind"]]), dimnames = list(NULL, v.attr_names, NULL))
  
  # starting values (first cycle)
  a.out[1,"TIME",]      <- 0
  a.out[1,"ALIVE",]     <- 1
  a.out[1,"AGE",]       <- round(qnorm(p = a.random[1,"AGE",], mean = l.inputs[["AGE_start_mean"]], sd = l.inputs[["AGE_start_sd"]]),0)
  a.out[1,"AGE",][a.out[1,"AGE",]<50] <- 50
  a.out[1,"AGE",][a.out[1,"AGE",]>99] <- 99
  a.out[1,"SEX",]       <- f.qcat(p_rand = a.random[1,"SEX",], p_cat = c(l.inputs[["p.SEX_start_male"]], 
                                                                         l.inputs[["p.SEX_start_female"]]), values = v.SEX_val)
  a.out[1,"EDU",]       <- f.qcat(p_rand = a.random[1,"EDU",], p_cat = l.inputs[["p.EDU_start"]], values = v.EDU_val)
  a.out[1,"RACEETH",]   <- f.qcat(p_rand = a.random[1,"RACEETH",], p_cat = l.inputs[["p.RACEETH_start"]], values = v.RACEETH_val)
  a.out[1,"INCOME",]    <- f.qcat(p_rand = a.random[1,"INCOME",], p_cat = l.inputs[["p.INCOME_start"]], values = v.INCOME_val)
  a.out[1,"MEDBUR",]    <- qunif(p = a.random[1,"APOE4",], min = min(v.MEDBUR_val), max = max(v.MEDBUR_val))
  a.out[1,"APOE4",]     <- f.qcat(p_rand = a.random[1,"APOE4",], p_cat = l.inputs[["p.APOE4_start"]], values = v.APOE4_val)
  
  
  a.out[1,"TX",]        <- 0
  # as.numeric(l.inputs[["Tx"]]==1)
  a.out[1,"SYN",]        <- 0    # Everyone starts healthy
  a.out[1,"BHA",]        <- NA   # Will be assigned as people are tested
  a.out[1,"CDR_track",]  <- NA   # Will be determined based on attributes
  a.out[1,"CDRfast_sd1",] <- qnorm(p = a.random[1,"CDRfast_sd1",], mean = 0, sd = l.inputs[["r.CDRfast_sd1"]])  # Random number
  a.out[1,"CDRslow_sd1",] <- qnorm(p = a.random[1,"CDRslow_sd1",], mean = 0, sd = l.inputs[["r.CDRslow_sd1"]])  # Random number
  a.out[1,"CDR",]       <- 0    # CDR-SB true score -- everyone starts cognitively normal
  a.out[1,"CDR_obs",]   <- NA   # CDR-SB observed score -- will be assigned as people are tested
  a.out[1,"SEV",]       <- NA   # Will be assigned as people are tested
  a.out[1,"SEV_obs",]   <- NA   # Will be assigned as people are tested
  a.out[1,"FUN",]       <- 0    # FAQ -- everyone starts with no functional impairment
  a.out[1,"BEH",]       <- NA
  a.out[1,"INSTIT",]    <- 0
  a.out[1,"QALY",]      <- NA
  a.out[1,"COST_care",] <- NA
  a.out[1,"COST_tx",]   <- NA
  
  # update progress bar
  setTxtProgressBar(pb, 1)
  
  # run subsequent cycles
  for(t in 2:l.inputs[["n.cycle"]]) {
    
    # TIME
    a.out[t,"TIME",] <- f.update_TIME(
      v.TIME.lag   = a.out[t-1,"TIME",]
    )
    
    alive.lag <- a.out[t-1,"ALIVE",]==1
    
    # ALIVE
    a.out[t,"ALIVE",alive.lag] <- f.update_ALIVE(
      v.AGE.lag     = a.out[t-1,"AGE",alive.lag], 
      v.SYN.lag     = a.out[t-1,"SYN",alive.lag],
      v.SEV.lag     = a.out[t-1,"SEV",alive.lag], 
      random_cycle  = a.random[t,"ALIVE",alive.lag], 
      m.lifetable   = l.inputs[["m.lifetable"]], 
      hr.mort_mci   = l.inputs[["hr.mort_mci"]], 
      hr.mort_mil   = l.inputs[["hr.mort_mil"]], 
      hr.mort_mod   = l.inputs[["hr.mort_mod"]], 
      hr.mort_sev   = l.inputs[["hr.mort_sev"]]
    )
    
    a.out[t,"ALIVE",!(alive.lag)] <- 0
    
    # identify those alive at current observation (to be used for subsetting the other functions 
    #       so they don't have to process the data of the individuals no longer alive)
    alive <- a.out[t,"ALIVE",]==1
    
    # number of individuals alive
    n.alive <- sum(alive)
    
    # AGE
    a.out[t,"AGE",alive] <- f.update_AGE(
      v.AGE.lag    = a.out[t-1,"AGE",alive]
    )
    
    # SEX
    a.out[t,"SEX",alive] <- f.update_SEX(
      v.SEX.lag = a.out[t-1,"SEX",alive]
    )
    
    # EDU
    a.out[t,"EDU",alive] <- f.update_EDU(
      v.EDU.lag = a.out[t-1,"EDU",alive]
    )
    
    # RACEETH
    a.out[t,"RACEETH",alive] <- f.update_RACEETH(
      v.RACEETH.lag = a.out[t-1,"RACEETH",alive]
    )
    
    # INCOME
    a.out[t,"INCOME",alive] <- f.update_INCOME(
      v.INCOME.lag = a.out[t-1,"INCOME",alive]
    )
    
    # MEDBUR
    a.out[t,"MEDBUR",alive] <- f.update_MEDBUR(
      v.MEDBUR.lag = a.out[t-1,"MEDBUR",alive]
    )
    
    # APOE4
    a.out[t,"APOE4",alive] <- f.update_APOE4(
      v.APOE4.lag = a.out[t-1,"APOE4",alive]
    )
    
    # MCI
    a.out[t,"SYN",alive] <- f.update_SYN(
      l.inputs        = l.inputs,
      v.AGE.lag       = a.out[t-1,"AGE",alive],
      v.EDU.lag       = a.out[t-1,"EDU",alive],
      v.APOE4.lag     = a.out[t-1,"APOE4",alive],
      v.MEDBUR.lag      = a.out[t-1,"MEDBUR",alive],
      v.INCOME.lag    = a.out[t-1,"INCOME",alive],
      v.SYN.lag       = a.out[t-1,"SYN",alive], 
      random_cycle    = a.random[t,"SYN",alive], 
      n.alive         = n.alive
    )
    
    # BHA
    a.out[t,"BHA",alive] <- f.update_BHA(
      v.SYN            = a.out[t,"SYN",alive],
      v.SYN.lag        = a.out[t-1,"SYN",alive],
      sens_BHA         = l.inputs[["sens_BHA"]],
      spec_BHA         = l.inputs[["spec_BHA"]],
      random_cycle     = a.random[t,"BHA",alive],
      n.alive          = n.alive
      )
    
    # CDR_track
    a.out[t,"CDR_track",alive] <- f.update_CDR_track(
      v.MEDBUR         = a.out[t,"MEDBUR",alive],
      random_cycle     = a.random[t,"CDR_track",alive],
      n.alive          = n.alive
    )
    
    # CDRfast_sd1, CDRslow_sd1
    a.out[t,"CDRfast_sd1",alive] <- a.out[t-1,"CDRfast_sd1",alive]
    a.out[t,"CDRslow_sd1",alive] <- a.out[t-1,"CDRslow_sd1",alive]
    
    # CDR (true)
    a.out[t,"CDR",alive] <- f.update_CDR(
      v.SYN            = a.out[t,"SYN",alive],
      v.SYN.lag        = a.out[t-1,"SYN",alive],
      cutoff_CDR       = l.inputs[["cutoff_CDR"]],
      v.CDR.lag   = a.out[t-1,"CDR",alive], 
      r.CDRfast_mean   = l.inputs[["r.CDRfast_mean"]],
      r.CDRslow_mean   = l.inputs[["r.CDRslow_mean"]],
      v.CDR_track      = a.out[t,"CDR_track",alive],
      v.CDRfast_sd1    = a.out[t,"CDRfast_sd1",alive],
      v.CDRslow_sd1    = a.out[t,"CDRslow_sd1",alive],
      r.CDR_sd2        = l.inputs[["r.CDR_sd2"]],
      random_cycle     = a.random[t,"CDR",alive],
      n.alive          = n.alive
    )
    
    
    # CDR_obs
    a.out[t,"CDR_obs",alive] <- f.update_CDR_obs(
      v.BHA            = a.out[t,"BHA",alive],
      v.CDR       = a.out[t,"CDR",alive],
      r.CDR_sd3        = l.inputs[["r.CDR_sd3"]],
      random_cycle     = a.random[t,"CDR_obs",alive],
      n.alive          = n.alive
    )

    # SEV (true)
    a.out[t,"SEV",alive] <- f.update_SEV(
      v.SYN          = a.out[t,"SYN",alive],
      v.CDR     = a.out[t,"CDR",alive],
      cutoff_CDR     = l.inputs[["cutoff_CDR"]],
      n.alive        = n.alive
    )
    
    
    # SEV_obs
    a.out[t,"SEV_obs",alive] <- f.update_SEV_obs(
      v.CDR_obs     = a.out[t,"CDR_obs",alive],
      cutoff_CDR     = l.inputs[["cutoff_CDR"]],
      n.alive        = n.alive
    ) 
    
    
    # FAQ 
    
    
    
    # TX
    a.out[t,"TX",alive] <- f.update_TX(
      n.alive              = n.alive, 
      v.TX.lag             = a.out[t-1,"TX",alive], 
      v.time               = a.out[t,"TIME",alive], 
      v.SYN                = a.out[t,"SYN",alive], 
      random_cycle         = a.random[t,"TX",alive], 
      Tx_t_max             = l.inputs[["Tx_t_max"]], 
      p.Tx_discontinuation = l.inputs[["p.Tx_discontinuation"]], 
      Tx                   = l.inputs[["Tx"]]
    )
    
    # update progress bar
    setTxtProgressBar(pb, t)
    
  }
  
  
  # The "round(score*2)/2" ensures that scores are increments of 0.5 
  a.out[,"CDR",] <- round(a.out[,"CDR",] * 2) / 2
  a.out[,"CDR_obs",] <- round(a.out[,"CDR_obs",] * 2) / 2
  
  # run time
  print(proc.time() - ptm)
  print(Sys.time() - stime)
  
  # return outcome
  return(a.out)
  
}


# apply QALYs and costs
f.qaly_cost <- function(a.out, l.inputs) {
  
  # QALY
  QALY0 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==0) * l.inputs[["u.healthy"]]
  QALY1 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==0) * l.inputs[["u.mci"]]
  QALY2 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==1) * l.inputs[["u.mil"]] +
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==2) * l.inputs[["u.mod"]] +
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==3) * l.inputs[["u.sev"]]
  QALY3 <- as.numeric(a.out[,"TIME",]==0 & a.out[,"TX",])
  QALY0[is.na(QALY0)] <- 0
  QALY1[is.na(QALY1)] <- 0
  QALY2[is.na(QALY2)] <- 0
  QALY3[is.na(QALY3)] <- 0
  
  # COST: treatment
  COST_tx <- as.numeric(a.out[,"TX",]) * l.inputs[["c.Tx"]]
  COST_tx[is.na(COST_tx)] <- 0
  
  # COST: care
  COST_care0 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==0) * l.inputs[["c.healthy"]]
  COST_care1 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==0) * l.inputs[["c.mci"]]
  COST_care2 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==1) * l.inputs[["c.mil"]] + 
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==2) * l.inputs[["c.mod"]] + 
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1 & a.out[,"SEV",]==3) * l.inputs[["c.sev"  ]]
  COST_care0[is.na(COST_care0)] <- 0
  COST_care1[is.na(COST_care1)] <- 0
  COST_care2[is.na(COST_care2)] <- 0
  
  # store
  a.out[,"QALY",] <- QALY0 + QALY1 + QALY2 + QALY3
  # a.out[,"COST_test",] <- NA
  a.out[,"COST_tx",] <- COST_tx
  # a.out[,"COST_fu",] <- NA
  a.out[,"COST_care",] <- COST_care0 + COST_care1 + COST_care2
  
  # return
  return(a.out)
  
}


# aggregate outcomes
f.out_aggregate <- function(a.out, l.inputs) {
  l.out <- vector(mode = "list", length = 0)
  
  # mean time alive
  l.out[["alive.trc"]] <- as.matrix(apply(X = a.out[,"ALIVE",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["alive.sum"]] <- sum(l.out[["alive.trc"]])
  l.out[["alive.trc.dis"]] <- as.matrix(f.discount(x = l.out[["alive.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["alive.sum.dis"]] <- sum(l.out[["alive.trc.dis"]])
  
  
  # mean time healthy
  l.out[["healthy.trc"]] <- as.matrix(apply(X = a.out[,"SYN",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["healthy.sum"]] <- sum(l.out[["healthy.trc"]])
  l.out[["healthy.trc.dis"]] <- as.matrix(f.discount(x = l.out[["healthy.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["healthy.sum.dis"]] <- sum(l.out[["healthy.trc.dis"]])
  l.out[["healthy_obs.trc"]] <- as.matrix(apply(X = a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  
  # time in MCI
  l.out[["MCI.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["MCI.sum"]] <- sum(l.out[["MCI.trc"]])
  l.out[["MCI.trc.dis"]] <- as.matrix(f.discount(x = l.out[["MCI.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["MCI.sum.dis"]] <- sum(l.out[["MCI.trc.dis"]])
  l.out[["MCI_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  
  # time in mild dementia
  l.out[["SEV1.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV1.sum"]] <- sum(l.out[["SEV1.trc"]])
  l.out[["SEV1.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV1.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV1.sum.dis"]] <- sum(l.out[["SEV1.trc.dis"]])
  l.out[["SEV1_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  
  # time in moderate dementia
  l.out[["SEV2.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV2.sum"]] <- sum(l.out[["SEV2.trc"]])
  l.out[["SEV2.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV2.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV2.sum.dis"]] <- sum(l.out[["SEV2.trc.dis"]])
  l.out[["SEV2_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  
  # time in severe dementia
  l.out[["SEV3.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV3.sum"]] <- sum(l.out[["SEV3.trc"]])
  l.out[["SEV3.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV3.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV3.sum.dis"]] <- sum(l.out[["SEV3.trc.dis"]])
  l.out[["SEV3_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  
  # temporary to select outcomes
  l.out[["mean_time_alive"]] <- l.out[["alive.sum.dis"]]
  l.out[["mean_time_healthy"]] <- l.out[["healthy.sum.dis"]]
  l.out[["mean_time_MCI"]] <- l.out[["MCI.sum.dis"]]
  l.out[["mean_time_SEV1"]] <- l.out[["SEV1.sum.dis"]]
  l.out[["mean_time_SEV2"]] <- l.out[["SEV2.sum.dis"]]
  l.out[["mean_time_SEV3"]] <- l.out[["SEV3.sum.dis"]]
  
  # time on treatment
  l.out[["mean_time_Tx"]] <- sum(a.out[,"TX",]==1, na.rm=TRUE)/l.inputs[["n.ind"]]
  
  # time in full-time care
  # l.out[["mean_time_FTC"]] <- sum(a.out[,"INSTIT",]==1, na.rm=TRUE)/l.inputs[["n.ind"]]
  
  # state trace (undiscounted) (true states)
  l.out[["state_trace"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 6, dimnames = list(NULL,c("healthy","mci","mil","mod","sev","dth")))
  l.out[["state_trace"]][,"healthy"] <- l.out[["healthy.trc"]]
  l.out[["state_trace"]][,"mci"] <- l.out[["MCI.trc"]]
  l.out[["state_trace"]][,"mil"] <- l.out[["SEV1.trc"]]
  l.out[["state_trace"]][,"mod"] <- l.out[["SEV2.trc"]]
  l.out[["state_trace"]][,"sev"] <- l.out[["SEV3.trc"]]
  l.out[["state_trace"]][,"dth"] <- apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]]
  # check rowsum
  rowSums(l.out[["state_trace"]])
  # trace institutionalized
  # l.out[["state_trace_instit"]] <- as.matrix((apply(X = a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]]))
  
  # state trace (undiscounted) (observed states)
  l.out[["state_trace_obs"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 6, dimnames = list(NULL,c("healthy","mci","mil","mod","sev","dth")))
  l.out[["state_trace_obs"]][,"healthy"] <- l.out[["healthy.trc"]]
  l.out[["state_trace_obs"]][,"mci"] <- l.out[["MCI_obs.trc"]]
  l.out[["state_trace_obs"]][,"mil"] <- l.out[["SEV1_obs.trc"]]
  l.out[["state_trace_obs"]][,"mod"] <- l.out[["SEV2_obs.trc"]]
  l.out[["state_trace_obs"]][,"sev"] <- l.out[["SEV3_obs.trc"]]
  l.out[["state_trace_obs"]][,"dth"] <- apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]]
  # check rowsum
  rowSums(l.out[["state_trace"]])
  
  
  # QALY
  l.out[["QALY"]] <- as.matrix(apply(X = a.out[,"QALY",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["QALY.sum"]] <- sum(l.out[["QALY"]])
  l.out[["QALY.dis"]] <- as.matrix(f.discount(x = l.out[["QALY"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["QALY.dis.sum"]] <- sum(l.out[["QALY.dis"]])
  
  # COST_test
  # l.out[["COST_test"]] <- as.matrix(apply(X = a.out[,"COST_test",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  # l.out[["COST_test.sum"]] <- sum(l.out[["COST_test"]])
  # l.out[["COST_test.dis"]] <- as.matrix(f.discount(x = l.out[["COST_test"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  # l.out[["COST_test.dis.sum"]] <- sum(l.out[["COST_test.dis"]])
  
  # COST_tx
  l.out[["COST_tx"]] <- as.matrix(apply(X = a.out[,"COST_tx",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["COST_tx.sum"]] <- sum(l.out[["COST_tx"]])
  l.out[["COST_tx.dis"]] <- as.matrix(f.discount(x = l.out[["COST_tx"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_tx.dis.sum"]] <- sum(l.out[["COST_tx.dis"]])
  
  # COST_fu
  # l.out[["COST_fu"]] <- as.matrix(apply(X = a.out[,"COST_fu",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  # l.out[["COST_fu.sum"]] <- sum(l.out[["COST_fu"]])
  # l.out[["COST_fu.dis"]] <- as.matrix(f.discount(x = l.out[["COST_fu"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  # l.out[["COST_fu.dis.sum"]] <- sum(l.out[["COST_fu.dis"]])
  # 
  # COST_care
  l.out[["COST_care"]] <- as.matrix(apply(X = a.out[,"COST_care",], MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["COST_care.sum"]] <- sum(l.out[["COST_care"]])
  l.out[["COST_care.dis"]] <- as.matrix(f.discount(x = l.out[["COST_care"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_care.dis.sum"]] <- sum(l.out[["COST_care.dis"]])
  
  # COST_tot
  l.out[["COST_tot"]]         <- l.out[["COST_tx"]]         + l.out[["COST_care"]]         
  # + l.out[["COST_test"]]         + l.out[["COST_fu"]]
  l.out[["COST_tot.sum"]]     <- l.out[["COST_tx.sum"]]     + l.out[["COST_care.sum"]]     
  # + l.out[["COST_test.sum"]]     + l.out[["COST_fu.sum"]]
  l.out[["COST_tot.dis"]]     <- l.out[["COST_tx.dis"]]     + l.out[["COST_care.dis"]]     
  # + l.out[["COST_test.dis"]]     + l.out[["COST_fu.dis"]]
  l.out[["COST_tot.dis.sum"]] <- l.out[["COST_tx.dis.sum"]] + l.out[["COST_care.dis.sum"]] 
  # + l.out[["COST_test.dis.sum"]] + l.out[["COST_fu.dis.sum"]]
  
  # net health benefit (NHB)
  l.out[["NHB"]]     <- l.out[["QALY.sum"]]     - l.out[["COST_tot.sum"]]     / 20000
  l.out[["NHB.dis"]] <- l.out[["QALY.dis.sum"]] - l.out[["COST_tot.dis.sum"]] / 20000
  
  # net monetary benefit (NMB)
  l.out[["NMB"]]     <- l.out[["QALY.sum"]]     * 20000 - l.out[["COST_tot.sum"]]
  l.out[["NMB.dis"]] <- l.out[["QALY.dis.sum"]] * 20000 - l.out[["COST_tot.dis.sum"]]
  
  # reporting table summed results
  l.out[["table_sum"]] <- matrix(
    data = c(
      l.out[["mean_time_alive"]],
      l.out[["mean_time_healthy"]],
      l.out[["mean_time_MCI"]],
      l.out[["mean_time_SEV1"]],
      l.out[["mean_time_SEV2"]],
      l.out[["mean_time_SEV3"]],
      l.out[["mean_time_Tx"]],
      # l.out[["mean_time_FTC"]],
      l.out[["QALY.sum"]],
      # l.out[["COST_test.sum"]],
      l.out[["COST_tx.sum"]],
      # l.out[["COST_fu.sum"]],
      l.out[["COST_care.sum"]],
      l.out[["COST_tot.sum"]],
      l.out[["NHB"]],
      l.out[["NMB"]],
      l.out[["QALY.dis.sum"]],
      # l.out[["COST_test.dis.sum"]],
      l.out[["COST_tx.dis.sum"]],
      # l.out[["COST_fu.dis.sum"]],
      l.out[["COST_care.dis.sum"]],
      l.out[["COST_tot.dis.sum"]],
      l.out[["NHB.dis"]],
      l.out[["NMB.dis"]]
    ),
    nrow = 1,
    ncol = 24-5,
    dimnames = list(NULL, c("mean_time_alive","mean_time_healthy","mean_time_MCI","mean_time_SEV1","mean_time_SEV2","mean_time_SEV3","mean_time_Tx",
                            # "mean_time_FTC",
                            "QALY.sum",
                            # "COST_test.sum",
                            "COST_tx.sum",
                            # "COST_fu.sum",
                            "COST_care.sum","COST_tot.sum","NHB","NMB",
                            "QALY.dis.sum",
                            # "COST_test.dis.sum",
                            "COST_tx.dis.sum",
                            # "COST_fu.dis.sum",
                            "COST_care.dis.sum","COST_tot.dis.sum","NHB.dis","NMB.dis")))
  
  # reporting table trace results
  l.out[["table_trace"]] <- matrix(
    data = c(
      l.out[["state_trace"]][,"healthy"],
      l.out[["state_trace"]][,"mci"],
      l.out[["state_trace"]][,"mil"],
      l.out[["state_trace"]][,"mod"],
      l.out[["state_trace"]][,"sev"],
      l.out[["state_trace"]][,"dth"],
      # l.out[["state_trace_instit"]],
      l.out[["QALY"]],
      # l.out[["COST_test"]],
      l.out[["COST_tx"]],
      # l.out[["COST_fu"]],
      l.out[["COST_care"]],
      l.out[["COST_tot"]],
      l.out[["QALY.dis"]],
      # l.out[["COST_test.dis"]],
      l.out[["COST_tx.dis"]],
      # l.out[["COST_fu.dis"]],
      l.out[["COST_care.dis"]],
      l.out[["COST_tot.dis"]]
    ),
    ncol = 19-5,
    dimnames = list(NULL, c("healthy","mci","mil","mod","sev","dth",
                            # "state_trace_instit",
                            "QALY",
                            # "COST_test",
                            "COST_tx",
                            # "COST_fu",
                            "COST_care","COST_tot","QALY.dis",
                            # "COST_test.dis",
                            "COST_tx.dis",
                            # "COST_fu.dis",
                            "COST_care.dis","COST_tot.dis"))
  )
  
  # return
  return(l.out)
  
}

# run strategies and incremental outcomes
f.out_summary <- function(a.random, l.inputs) {
  
  # random values
  if(is.null(a.random)) {
    a.random <- array(data = NA, dim = c(l.inputs[["n.cycle"]], n.attr, l.inputs[["n.ind"]]), dimnames = list(NULL, v.attr_names, NULL))
    a.random[,,] <- runif(n = length(a.random))
  }
  
  # strategy 1
  l.inputs_strat1 <- l.inputs
  l.inputs_strat1[["strategy"]] <- l.inputs[["strategy_strat1"]]
  l.inputs_strat1[["Tx"]] <- l.inputs[["Tx_strat1"]]
  a.out_strat1 <- f.run(a.random = a.random, l.inputs = l.inputs_strat1)
  a.out_qc_strat1 <- f.qaly_cost(a.out = a.out_strat1, l.inputs = l.inputs_strat1)
  out_strat1 <- f.out_aggregate(a.out = a.out_qc_strat1, l.inputs = l.inputs_strat1)
  
  # strategy 2
  l.inputs_strat2 <- l.inputs
  l.inputs_strat2[["strategy"]] <- l.inputs[["strategy_strat2"]]
  l.inputs_strat2[["Tx"]] <- l.inputs[["Tx_strat2"]]
  a.out_strat2 <- f.run(l.inputs = l.inputs_strat2, a.random = a.random)
  a.out_qc_strat2 <- f.qaly_cost(a.out = a.out_strat2, l.inputs = l.inputs_strat2)
  out_strat2 <- f.out_aggregate(a.out = a.out_qc_strat2, l.inputs = l.inputs_strat2)
  
  # summary outcomes
  m.out <- matrix(
    data = NA,
    nrow = 5,
    ncol = ncol(out_strat1[["table_sum"]]),
    dimnames = list( c("strategy 1 (cau)","strategy 2 (dmt)","strategy 3 (symptomatic)","incr. strategy 2-1","incr. strategy 3-1"), colnames(out_strat1[["table_sum"]]) )
  )
  m.out[1,] <- out_strat1[["table_sum"]]
  m.out[2,] <- out_strat2[["table_sum"]]
  
  m.out[4,] <- m.out[2,] - m.out[1,]
  m.out[5,] <- m.out[3,] - m.out[1,]
  m.out <- cbind(m.out, ICER = NA, ICER.dis = NA)
  m.out[4,"ICER"] <- (m.out[2,"COST_tot.sum"    ] - m.out[1,"COST_tot.sum"    ]) / (m.out[2,"QALY.sum"    ] - m.out[1,"QALY.sum"    ])
  m.out[4,"ICER.dis"] <- (m.out[2,"COST_tot.dis.sum"] - m.out[1,"COST_tot.dis.sum"]) / (m.out[2,"QALY.dis.sum"] - m.out[1,"QALY.dis.sum"])
  
  return(list(
    a.random = a.random,
    l.inputs = l.inputs,
    a.out_strat1 = a.out_strat1, a.out_strat2 = a.out_strat2,
    a.out_qc_strat1 = a.out_qc_strat1, a.out_qc_strat2 = a.out_qc_strat2,
    out_strat1 = out_strat1, out_strat2 = out_strat2,
    m.out = m.out
  ))
  
}

f.figures <- function(l.out, l.inputs) {
  
  # Reshape syndrome data to long format for use with ggplot.
  t.trace_syndrome <- as.data.frame(l.out$table_trace[ ,c("healthy", "mci", "mil", "mod", "sev", "dth")],) %>%
    mutate(year = 1:l.inputs[["n.cycle"]])
  t.trace_syndrome_long <- t.trace_syndrome %>%
    pivot_longer(cols = -year, names_to = "syndrome", values_to = "proportion") %>%
    mutate(syndrome = fct_rev(factor(syndrome, levels = c("healthy", "mci", "mil", "mod", "sev", "dth"))))
  
  fig.syndrome_trace <- ggplot(t.trace_syndrome_long, aes(x = year, y = proportion, fill = syndrome)) +
    geom_area(alpha = 0.8, position = "stack") + 
    geom_path(aes(group = syndrome), position = "stack", color = "black", linewidth = 0.5) + 
    scale_x_continuous(breaks = unique(t.trace_syndrome_long$year), minor_breaks = NULL) +
    scale_fill_manual(values = c("white","orange","purple","green","yellow","pink"),
                      labels = c("Death","Severe dementia","Moderate dementia","Mild dementia","MCI","Cognitively intact")) +
    labs(title = "Progression of Cognitive Impairment",
         x = "Year",
         y = "Proportion of Population",
         fill = "Cognitive Status")
  
  return(fig.syndrome_trace)
}

######################################## 2.4. FUNCTIONS: WRAPPERS ########################################

f.wrap_run <- function(a.random, l.inputs) {
  output <- f.run(a.random = a.random, l.inputs = l.inputs)
  aggregated_results <- f.out_aggregate(a.out = output, l.inputs = l.inputs)
  figures <- f.figures(l.out = aggregated_results, l.inputs = l.inputs)
  
  print(figures)
  
  invisible(return(list(
    output = output,
    aggregated_results = aggregated_results,
    figures = figures)
  ))
}


######################################## 2.5. RUN MODEL ########################################

model <- f.wrap_run(a.random = a.random, l.inputs = l.inputs)

