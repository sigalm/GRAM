# ****************************************************************
# ======= 2. GRAM HELPER FUNCTIONS ======= 
# ****************************************************************

######################################## 2.1. HELPER FUNCTIONS: GENERIC ########################################
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
    if (any(p_cat == 0) | any(p_cat == 1)) {        # if p_cat contains absolute probabilities, offset them to avoid errors
      p_cat[p_cat == 0] <- 1e-10
      p_cat[p_cat == 1] <- 1 - 1e-10
    }
    
    # function
    breaks <- c(0,cumsum(p_cat)) 
    breaks[duplicated(breaks)] <- breaks[duplicated(breaks)] + 1e-10 
    out <- cut(x = p_rand, breaks = breaks, labels = FALSE, include.lowest = TRUE) # categorize random values based on user-defined breaks
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

######################################## discounting

f.discount <- function(x, l.inputs, discount_rate, n.cycle) {
  # x can be scalar, vector or matrix
  as.matrix(x) / (1 + discount_rate)^(0:(n.cycle - 1))
}


######################################## 2.2. HELPER FUNCTIONS: ATTRIBUTE UPDATE ########################################
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
    alive.lag,
    v.AGE.lag, v.SYN.lag, v.SEV.lag, 
    m.lifetable, 
    hr.mort_mci, hr.mort_mil, hr.mort_mod, hr.mort_sev, 
    hr.mort_age,
    random_cycle) 
{
  
  # generate life table coordinates for looking up age-specific mortality 
  # (see https://adv-r.hadley.nz/subsetting.html paragraph 4.2.3 subsetting > selecting multiple elements > subsetting)
  lifetable_lookup_coordinates <- matrix(data = round(v.AGE.lag,0), ncol = 1) - 50 + 1
  
  # determine relative mortality risk related to syndrome and severity
  healthy <- v.SYN.lag==0
  mci <- v.SYN.lag==1 & v.SEV.lag==0
  mil <- v.SYN.lag==1 & v.SEV.lag==1
  mod <- v.SYN.lag==1 & v.SEV.lag==2
  sev <- v.SYN.lag==1 & v.SEV.lag==3
  
  v.risk_ratio <- rep(NA, length(v.AGE.lag))
  v.risk_ratio[healthy] <- 1
  v.risk_ratio[mci] <- hr.mort_mci
  v.risk_ratio[mil] <- hr.mort_mil 
  v.risk_ratio[mod] <- hr.mort_mod * ((v.AGE.lag[mod]<70) + hr.mort_age[2] * (v.AGE.lag[mod]>=70 & v.AGE.lag[mod]<80) + hr.mort_age[3] * (v.AGE.lag[mod]>=80))
  v.risk_ratio[sev] <- hr.mort_sev * ((v.AGE.lag[sev]<70) + hr.mort_age[2] * (v.AGE.lag[sev]>=70 & v.AGE.lag[sev]<80) + hr.mort_age[3] * (v.AGE.lag[sev]>=80))
  
  
  # calculate death probability using life table, relative risk and adjustment for cycle length
  prob_death <- pmin(f.adjustprobability(
    p = m.lifetable[lifetable_lookup_coordinates], 
    t_new = 1, 
    t_old = 1,
    RR = v.risk_ratio),1)
  
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
  
  prob_none_to_mci <- (1 - exp(-hazard)) * l.inputs[["rr.Px_mci"]]
  
  # apply dependent on previous state
  symptoms[v.SYN.lag == 0] <- as.numeric(prob_none_to_mci[v.SYN.lag==0] > random_cycle[v.SYN.lag==0])
  symptoms[v.SYN.lag == 1] <- 1
  
  return(symptoms)
}

######################################## CDR-SB fast/slow track

f.update_CDR_track <- function(v.SEV.lag, n.alive) {
  
  track <- rep(NA, n.alive)
  
  # People with true MCI will be slow
  # People with true dementia will be fast
  
  track[v.SEV.lag == 0] <- 0
  track[v.SEV.lag >= 1] <- 1 
  
  return(track)
}

######################################## CDR-SB - true disease status

f.update_CDR <- function(v.SYN, v.SYN.lag, cutoff_CDR, v.CDR.lag, 
                         r.CDRfast_mean, r.CDRslow_mean, v.CDR_track,
                         v.CDRfast_sd1, v.CDRslow_sd1, r.CDR_sd2, 
                         v.TX.lag, rr.Tx_mci, random_cycle, n.alive) {
  
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
  
  delta <- 
    r.CDRfast_mean * v.CDR_track + r.CDRslow_mean * (1-v.CDR_track) +         # Add mean increase in CDR-SB score
    v.CDRfast_sd1 * v.CDR_track + v.CDRslow_sd1 * (1-v.CDR_track) +           # Add individual-level deviation from mean
    qnorm(p = random_cycle, mean = 0, sd = r.CDR_sd2)                         # Add within-individual deviation
  
  delta_tx <- delta * v.TX.lag * rr.Tx_mci + delta * !v.TX.lag * 1
  
  cdr[mci_still] <- pmin(v.CDR.lag[mci_still] + delta_tx[mci_still],
                         cutoff_CDR["max"])
  
  return(cdr)
}

######################################## SEV

f.update_SEV <- function(v.SYN, v.CDR, cutoff_CDR, n.alive) {
  
  # Assign dementia severity based on cognitive score
  sev <- rep(NA, n.alive)
  
  sev[v.SYN == 1 & v.CDR < cutoff_CDR["mci"]] <- 0
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["mci"], cutoff_CDR["mild"])] <- 0
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["mild"], cutoff_CDR["moderate"])] <- 1
  sev[v.SYN == 1 & between(v.CDR, cutoff_CDR["moderate"], cutoff_CDR["severe"])]<- 2
  sev[v.SYN == 1 & v.CDR >= cutoff_CDR["severe"]] <- 3
  
  return(sev)
}


######################################## BHA

f.update_BHA <- function(v.SYN, v.BHA.lag, v.SEV, sens_BHA, spec_BHA, random_cycle, n.alive) {
  
  bha <- rep(-9, n.alive)
  
  prior_pos <- v.BHA.lag == 1
  
  healthy <- !prior_pos & v.SYN == 0
  bha[healthy] <- as.numeric((1 - spec_BHA) > random_cycle[healthy])
  
  mci <- !prior_pos & v.SYN == 1 & v.SEV == 0
  bha[mci] <- as.numeric(sens_BHA[1] > random_cycle[mci])
  
  dementia <- !prior_pos & v.SYN == 1 & v.SEV >= 1
  bha[dementia] <- as.numeric(sens_BHA[2] > random_cycle[dementia])
  
  bha[prior_pos & v.SYN == 0] <- as.numeric((1 - spec_BHA) > random_cycle[prior_pos & v.SYN == 0])
  bha[prior_pos & v.SYN == 1] <- 1
  
  return(bha)
}


######################################## CDR-SB - observed disease status

f.update_CDR_obs <- function(v.BHA, v.CDR, r.CDR_sd3, random_cycle, n.alive) {
  
  cdr_obs <- rep(-9, n.alive)
  
  cdr_obs[v.BHA == 1] <- v.CDR[v.BHA == 1] + 
    qnorm(p = random_cycle[v.BHA == 1], mean = 0, sd = r.CDR_sd3)        # Add measurement error (rater reliability)
  
  return(cdr_obs)
}


######################################## SEV_obs

f.update_SEV_obs <- function(v.CDR_obs, cutoff_CDR, n.alive) {
  
  # Assign dementia severity based on cognitive score
  sev_obs <- rep(-9, n.alive)
  
  sev_obs[between(v.CDR_obs, cutoff_CDR["healthy"], cutoff_CDR["mci"])] <- 0
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

f.update_TX <- function(v.TX.lag, v.SEV_obs, v.SEV_obs.lag, Tx_counter, n.alive, random_cycle, Tx_t_max, p.Tx, Tx) {
  
  # start with empty vector
  tx <- rep(0, n.alive)
  
  # previous Tx state * previous syndrome==MCI * time shorter than maximum treatment duration
  if(Tx==1) {
    
    initiate_tx <- v.SEV_obs.lag == -9 & (v.SEV_obs == 0 | v.SEV_obs == 1)
    tx[initiate_tx] <- f.qcat(p_rand = random_cycle[initiate_tx], p_cat = p.Tx, values = c(0,1))
    
    continue_tx <- v.TX.lag == 1 & v.SEV_obs < 2
    tx[continue_tx] <- v.TX.lag[continue_tx] * as.numeric(Tx_counter[continue_tx] < Tx_t_max)
  }
  
  
  return(tx)
  
}


######################################## 2.3. HELPER FUNCTIONS: FORMAT OUTPUTS ########################################

format_reside_time_table <- function(reside_time_data, scenario_title = NULL) {
  table <- flextable(reside_time_data) %>%
    set_header_labels(
      age_group = "Age Group",
      mci = "MCI",
      mil = "Mild Dementia",
      mod = "Moderate Dementia",
      sev = "Severe Dementia"
    ) %>%
    add_header_row(
      values = c("", "Average Reside Time in Years, by Age at Onset"),
      colwidths = c(1, 4)
    ) %>%
    add_header_lines(values = scenario_title) %>%
    bold(j = 1) %>%
    theme_vanilla() %>%
    autofit() %>%
    align(j = 2:5, align = "center", part = "all") %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all")
  
  return(table) 
}

format_prevalence_table <- function(prevalence_data, denom = 1000, scenario_title = NULL) {
  
  dat <- prevalence_data
  dat[ ,-1] <- round(dat[ ,-1] * denom, digits = 2)
  
  table <- flextable(dat) %>%
    set_header_labels(
      age_group = "Age",
      avg_mci = "MCI",
      avg_mil = "Mild Dementia",
      avg_mod = "Moderate Dementia",
      avg_sev = "Severe Dementia"
    ) %>%
    add_header_row(
      values = c("", paste0("Prevalence by Severity by Age, per ", denom)),
      colwidths = c(1, 4)
    ) %>%
    add_header_lines(values = scenario_title) %>%
    bold(j = 1) %>%
    theme_vanilla() %>%
    autofit() %>%
    align(j = 2:5, align = "center", part = "all") %>%
    bold(part = "header") %>%
    fontsize(size = 10, part = "all")
  
}
