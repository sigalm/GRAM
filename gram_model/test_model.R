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
v.attr_names <- c("TIME","ALIVE","AGE","SEX","TX","SYN","COG","SEV", "FUN","EDU","RACEETH","OBES","AGE_MCI",
                  "BEH","INSTIT","QALY","COST_care", "COST_tx")
n.attr <- length(v.attr_names) # number        of attributes

# define or describe possible attribute values
v.ALIVE_val <- c(0,1)    # 0 = dead, 1 = alive
v.SEX_val <- c(1,2)      # 1 = male, 2 = female
v.EDU_val <- c(1,2,3)    # 1 = college or more, 2 = high school or GED, 3 = less than high school 
v.RACEETH_val <- c(0,1) # 0 = not Hispanic or Black, 1 = Hispanic or Black
v.OBES_val <- c(0,1)      # 0 = not obese, 1 = obese
v.APOE4_val <- c(0,1,2)    # 0 = non-carrier, 1 = hetero, 2 = homo
v.SYN_val <- c(0,1,2)     # 0 = healthy, 1 = MCI, 2 = dementia
v.SEV_val <- c(1,2,3)    # 1 = mild, 2 = moderate, 3 = severe
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
l.inputs[["seed_pa"]] <- 12345 # seed for generating random values that drive probabilistic analysis
l.inputs[["n.psa"]] <- 10 # number of PSA iterations


######################################## 1.2. EXTERNAL MODEL INPUTS ########################################
# model inputs for deterministic analysis
l.inputs[["AGE_start_mean"]] <- 50
l.inputs[["AGE_start_sd"]] <- 0

l.inputs[["p.SEX_start_male"]] <- 0.5
l.inputs[["p.SEX_start_female"]] <- 0.5

l.inputs[["p.EDU_start"]] <- c(0.4, 0.5, 0.1) # p for college, high school, less than high school, respectively. must add to 1.
l.inputs[["p.OBES_start"]] <- c(0.95, 0.05) # p for OBES = 0 and OBES = 1, respectively
l.inputs[["p.RACEETH_start"]] <- c(0.65, 0.35) # p for RACEETH = 0 (not Hisp/Black) and RACEETH = 1 (Hisp/Black)

l.inputs[["hr.mort_mci"]] <- 1.82
l.inputs[["hr.mort_mil"]] <- 2.92
l.inputs[["hr.mort_mod"]] <- 3.85
l.inputs[["hr.mort_sev"]] <- 9.52

l.inputs[["logit_intercept"]] <- -5.0
l.inputs[["logit_b1"]] <- 0.001
l.inputs[["logit_b11"]] <- -0.0001
l.inputs[["logit_b12"]] <- 0.00005
l.inputs[["logit_b2"]] <- 0.001
l.inputs[["logit_b3"]] <- 0.001
l.inputs[["logit_b4"]] <- 0.001
l.inputs[["logit_b1EDU"]] <- -0.002
l.inputs[["logit_b1RACEETH"]] <- -0.002
l.inputs[["logit_b1OBES"]] <- -0.002

l.inputs[["weibull_mu"]] <- -1.53141
l.inputs[["weibull_gamma"]] <- 0.757504

l.inputs[["p.land_dist"]] <- c(0.75,0.20,0.05)  # PLACEHOLER p for landing at mild, moderate, or severe dementia

# l.inputs[["instit_constant"]] <- -4.255301
# l.inputs[["instit_COG2"]] <- 0.2345492
# l.inputs[["instit_COG3"]] <- 0.3538726
# l.inputs[["instit_FUN2"]] <- 0.7485477
# l.inputs[["instit_FUN3"]] <- 1.4328740
# l.inputs[["instit_BEH2"]] <- 0.2494476
# l.inputs[["instit_BEH3"]] <- 0.6015229

l.inputs[["u.healthy"]] <- 0.85 # PLACEHOLDER - do we want to apply age-related decreases to this?
l.inputs[["u.mci"]] <- 0.73
l.inputs[["u.mil"]] <- 0.69
l.inputs[["u.mod"]] <- 0.53
l.inputs[["u.sev"]] <- 0.38

# l.inputs[["c.mci_home"]] <- 13364
# l.inputs[["c.mil_home"]] <- 26727
# l.inputs[["c.mod_home"]] <- 31644
# l.inputs[["c.sev_home"]] <- 40645
# l.inputs[["c.mil_instit"]] <- 111902
# l.inputs[["c.mod_instit"]] <- 111902
# l.inputs[["c.sev_instit"]] <- 113523
l.inputs[["c.healthy"]] <- 0 # PLACEHOLDER - do we want to apply age-related decreases to this?
l.inputs[["c.mci"]] <- 13364
l.inputs[["c.mil"]] <- 26727
l.inputs[["c.mod"]] <- 31644
l.inputs[["c.sev"]] <- 40645
l.inputs[["c.Tx"]] <- 5000

l.inputs[["rr.Tx_mci"]] <- 0.8
l.inputs[["Tx_t_max"]] <- 5
l.inputs[["r.discount_QALY"]] <- 0.03
l.inputs[["r.discount_COST"]] <- 0.03
l.inputs[["m.lifetable"]] <- array(data = read.csv("../gram_data/table_prob_die_next_year.csv")[ , 2], dim = c(51,1),
                                   dimnames = list(50:100, "q"))




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

######################################## OBES

f.update_OBES <- function(v.OBES.lag) {
  obes <- v.OBES.lag
  return(obes)
}

######################################## EDU

f.update_EDU <- function(v.EDU.lag) {
  edu <- v.EDU.lag
  return(edu)
}

######################################## SYN

f.update_SYN <- function(l.inputs, v.AGE.lag, v.EDU.lag, v.RACEETH.lag, v.OBES.lag, 
                         v.time, v.time.lag, mu, gamma, v.SYN.lag, v.TX.lag, rr.Tx_mci, random_cycle, n.alive) {
  
  # start with empty vector
  symptoms <- rep(NA, n.alive)
  
  # logistic regression (determine healthy --> MCI)
  x <- l.inputs[["logit_intercept"]] + 
    l.inputs[["logit_b1"]] * (v.AGE.lag - 50) + 
    l.inputs[["logit_b2"]] * v.EDU.lag + 
    l.inputs[["logit_b3"]] * v.RACEETH.lag + 
    l.inputs[["logit_b4"]] * v.OBES.lag + 
    l.inputs[["logit_b11"]] * (v.AGE.lag - 50)^2 + 
    l.inputs[["logit_b12"]] * (v.AGE.lag - 50)^3 + 
    l.inputs[["logit_b1EDU"]] * (v.AGE.lag - 50) * v.EDU.lag  + 
    l.inputs[["logit_b1RACEETH"]] * (v.AGE.lag - 50) * v.RACEETH.lag + 
    l.inputs[["logit_b1OBES"]] * (v.AGE.lag - 50) * v.OBES.lag
  
  TP_none_to_mci <- exp(x) / (1 + exp(x))    
  
  # survival function (determine MCI --> dementia)
  temp <- (1 - f.predict_weibull(time = v.time, mu = mu, gamma = gamma) / f.predict_weibull(time = v.time.lag, mu = mu, gamma = gamma))
  
  # RR related to treatment effect
  TP_mci_to_dementia <- temp * v.TX.lag * rr.Tx_mci + temp * !v.TX.lag * 1
  
  # apply dependent on previous state
  symptoms[v.SYN.lag==0] <- as.numeric(TP_none_to_mci[v.SYN.lag==0] > random_cycle[v.SYN.lag==0])
  symptoms[v.SYN.lag==1] <- as.numeric(TP_mci_to_dementia[v.SYN.lag==1] > random_cycle[v.SYN.lag==1]) * 2
  symptoms[v.SYN.lag==1 & symptoms!=2] <- 1
  symptoms[v.SYN.lag==2] <- 2
  
  # return
  return(symptoms)
}


######################################## CDR-SB

f.update_COG <- function(v.SYN, v.SYN.lag, v.COG.lag, p.land_dist, random_cycle, n.alive) {
  #discuss the issue of median vs. mean rate of change
  cog <- rep(NA, n.alive)
  
  healthy <- v.SYN==0
  cog[healthy] <- 0
  
  # Assign initial CDR-SB score for those who develop MCI in the current cycle
  # Cut-off values from: O'Bryant et al 2012 (PMC3409562) Table 2
  mci_new <- v.SYN==1 & v.SYN.lag==0
  cog[mci_new] <- runif(n=sum(mci_new), min = 0.5, max = 4.0)
  
  # Update CDR-SB score for those who remain in MCI
  mci_still <- v.SYN==1 & v.SYN.lag==1
  cog[mci_still]  <- pmin(v.COG.lag[mci_still] * (1 + 1.67/18), 4.0)
  
  # Assign CDR-SB score for those who develop dementia in the current cycle
  dementia_new <- v.SYN==2 & v.SYN.lag!=2
  mild <- dementia_new & (p.land_dist[1] > random_cycle)
  moderate <- dementia_new & between(random_cycle, p.land_dist[1], (p.land_dist[1] + p.land_dist[2]))
  severe <- dementia_new & (random_cycle > (p.land_dist[1] + p.land_dist[2]))
    
  cog[mild] <- runif(n=sum(mild), min = 4.5, max = 9.0)
  cog[moderate] <-runif(n=sum(moderate), min = 9.5, max = 15.5)
  cog[severe] <- runif(n=sum(severe), min = 16.0, max = 18.0)
  
  #Update CDR-SB score for those who remain in dementia
  dementia_still <- v.SYN==2 & v.SYN.lag==2
  cog[dementia_still] <- pmin(v.COG.lag[dementia_still] * (1 + 2/18), 18.0)
  
  # return
  return(cog)
}

######################################## SEV

f.update_SEV <- function(v.SYN, v.COG, n.alive) {
  
  # Assign dementia severity based on cognitive score
  dementia <- v.SYN==2
  sev <- rep(NA, n.alive)
  
  sev[v.COG >= 4.5 & v.COG < 9.5] <- 1
  sev[v.COG >= 9.5 & v.COG < 16.0] <- 2
  sev[v.COG >= 16.0] <- 3
  
  return(sev)
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
  a.out[1,"OBES",]      <- f.qcat(p_rand = a.random[1,"OBES",], p_cat = l.inputs[["p.OBES_start"]], values = v.OBES_val)
  a.out[1,"RACEETH",]  <- f.qcat(p_rand = a.random[1,"RACEETH",], p_cat = l.inputs[["p.RACEETH_start"]], values = v.RACEETH_val)
  
  a.out[1,"TX",]        <- 0
  # as.numeric(l.inputs[["Tx"]]==1)
  a.out[1,"SYN",]       <- 0
  a.out[1,"SEV",]       <- NA    # Will be assigned as people develop dementia
  a.out[1,"COG",]       <- 0    # CDR-SB -- everyone starts cognitively normal
  a.out[1,"FUN",]       <- 0     # FAQ -- everyone starts with no functional impairment
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
    
    # OBES
    a.out[t,"OBES",alive] <- f.update_OBES(
      v.OBES.lag = a.out[t-1,"OBES",alive]
    )
    
    # RACEETH
    a.out[t,"RACEETH",alive] <- f.update_RACEETH(
      v.RACEETH.lag = a.out[t-1,"RACEETH",alive]
    )
    
    # SYN
    a.out[t,"SYN",alive] <- f.update_SYN(
      l.inputs        = l.inputs,
      v.AGE.lag       = a.out[t-1,"AGE",alive],
      v.EDU.lag       = a.out[t-1,"EDU",alive],
      v.RACEETH.lag   = a.out[t-1,"RACEETH",alive],
      v.OBES.lag      = a.out[t-1,"OBES",alive],
      v.time.lag      = a.out[t-1,"TIME",alive], 
      v.SYN.lag       = a.out[t-1,"SYN",alive], 
      v.TX.lag        = a.out[t-1,"TX",alive], 
      v.time          = a.out[t,"TIME",alive], 
      random_cycle    = a.random[t,"SYN",alive], 
      mu              = l.inputs[["weibull_mu"]],
      gamma           = l.inputs[["weibull_gamma"]], 
      rr.Tx_mci       = l.inputs[["rr.Tx_mci"]],
      n.alive         = n.alive
    )
    
    # COG
    a.out[t,"COG",alive] <- f.update_COG(
      n.alive          = n.alive, 
      v.SYN            = a.out[t,"SYN",alive],
      v.SYN.lag        = a.out[t-1,"SYN",alive],
      v.COG.lag        = a.out[t-1,"COG",alive],
      p.land_dist      = l.inputs[["p.land_dist"]], 
      random_cycle     = a.random[t,"COG",alive]
    )
    
    # SEV
    a.out[t,"SEV",alive] <- f.update_SEV(
      v.SYN          = a.out[t,"SYN",alive],
      v.COG          = a.out[t,"COG",alive],
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
  a.out[,"COG",] <- round(a.out[,"COG",] * 2) / 2
  
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
  QALY1 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1) * l.inputs[["u.mci"]]
  QALY2 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==1) * l.inputs[["u.mil"]] +
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==2) * l.inputs[["u.mod"]] +
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==3) * l.inputs[["u.sev"]]
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
  COST_care1 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==1) * l.inputs[["c.mci"]]
  COST_care2 <- as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==1) * l.inputs[["c.mil"]] + 
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==2) * l.inputs[["c.mod"]] + 
    as.numeric(a.out[,"ALIVE",] & a.out[,"SYN",]==2 & a.out[,"SEV",]==0) * l.inputs[["c.sev"  ]]
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
  
  # time in MCI
  l.out[["MCI.trc"]] <- as.matrix(apply(X = a.out[,"SYN",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["MCI.sum"]] <- sum(l.out[["MCI.trc"]])
  l.out[["MCI.trc.dis"]] <- as.matrix(f.discount(x = l.out[["MCI.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["MCI.sum.dis"]] <- sum(l.out[["MCI.trc.dis"]])
  
  # time in mild dementia
  l.out[["SEV1.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV1.sum"]] <- sum(l.out[["SEV1.trc"]])
  l.out[["SEV1.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV1.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV1.sum.dis"]] <- sum(l.out[["SEV1.trc.dis"]])
  
  # time in moderate dementia
  l.out[["SEV2.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV2.sum"]] <- sum(l.out[["SEV2.trc"]])
  l.out[["SEV2.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV2.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV2.sum.dis"]] <- sum(l.out[["SEV2.trc.dis"]])
  
  # time in severe dementia
  l.out[["SEV3.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV3.sum"]] <- sum(l.out[["SEV3.trc"]])
  l.out[["SEV3.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV3.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV3.sum.dis"]] <- sum(l.out[["SEV3.trc.dis"]])
  
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
  
  # state trace (undiscounted)
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

