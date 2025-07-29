######################################## GRAM MODULE: HEALTH STATUS ########################################
# This section defines the functions for updating health status variables.

#### Module Wrapper ####
f.module_true_health <- function(a.out, t, a.random, alive, n.alive) {
  
  # MEDBUR
  a.out[t,"MEDBUR",alive] <- f.update_MEDBUR(
    v.MEDBUR.lag = a.out[t-1,"MEDBUR",alive],
    v.AGE.lag    = a.out[t-1, "AGE", alive],
    coef_MEDBUR  = l.inputs[["coef_MEDBUR"]],
    amplification = l.inputs[["amplification_MEDBUR"]],
    max_MEDBUR   = max(l.inputs[["v.MEDBUR_val"]]),
    random_cycle = a.random[t, "MEDBUR", alive],
    n.alive      = n.alive
  )
  
  # APOE4
  a.out[t,"APOE4",alive] <- f.update_APOE4(
    v.APOE4.lag = a.out[t-1,"APOE4",alive]
  )
  
  # SYN
  if (t <= l.inputs[["n.cycle"]] - 2) {
    a.out[t,"SYN",alive] <- f.update_SYN(
      l.inputs        = l.inputs,
      v.AGE.tplus2    = a.out[t,"AGE",alive] + 2,
      v.SEX.lag       = a.out[t-1,"SEX",alive],
      v.RACEETH.lag   = a.out[t-1,"RACEETH",alive],
      v.EDU.lag       = a.out[t-1,"EDU",alive],
      v.APOE4.lag     = a.out[t-1,"APOE4",alive],
      v.MEDBUR.lag    = a.out[t-1,"MEDBUR",alive],  
      v.INCOME.lag    = a.out[t-1,"INCOME",alive],
      v.SYN.lag       = a.out[t-1,"SYN",alive],
      v.SYN.lag2      = a.out[t-2,"SYN",alive],
      v.MEMLOSS.lag   = a.out[t-1,"MEMLOSS",alive],
      random_tplus2   = a.random[t+2,"SYN",alive], 
      n.alive         = n.alive
    )} else {
      a.out[t,"SYN",alive] <- a.out[t-1,"SYN",alive]    # last two cycles no change in syndrome            
    }
  
  # MEMLOSS 
  a.out[t,"MEMLOSS",alive] <- f.update_MEMLOSS(
    v.MEMLOSS.lag = a.out[t-1,"MEMLOSS",alive],
    v.SYN         = a.out[t,"SYN",alive],
    v.SYN.lag     = a.out[t-1,"SYN",alive],
    v.AGE         = a.out[t,"AGE",alive],
    v.EDU.lag     = a.out[t-1,"EDU",alive],
    v.SEX.lag     = a.out[t-1,"SEX",alive],
    v.RACEETH.lag = a.out[t-1,"RACEETH",alive],
    v.APOE4.lag   = a.out[t-1,"APOE4",alive],
    v.MEDBUR.lag  = a.out[t-1,"MEDBUR",alive],
    v.INCOME.lag  = a.out[t-1,"INCOME",alive],
    l.inputs      = l.inputs,
    p.MEMLOSS_new = l.inputs[["p.MEMLOSS_new"]],
    random_cycle  = a.random[t,"MEMLOSS",alive],
    n.alive       = n.alive 
  )
  
  # CDR_track
  a.out[t,"CDR_track",alive] <- f.update_CDR_track(
    v.SEV.lag        = a.out[t-1,"SEV",alive],
    n.alive          = n.alive
  )
  
  # CDRfast_sd1, CDRslow_sd1
  a.out[t,"CDRfast_sd1",alive] <- a.out[t-1,"CDRfast_sd1",alive]
  a.out[t,"CDRslow_sd1",alive] <- a.out[t-1,"CDRslow_sd1",alive]
  
  if (length(l.inputs[["r.CDRfast_mean"]]) > 1) {
    r.CDRfast_mean <- l.inputs[["r.CDRfast_mean"]][t]
  } else {
    r.CDRfast_mean <- l.inputs[["r.CDRfast_mean"]]
  }
  
  if (length(l.inputs[["r.CDRslow_mean"]]) > 1) {
    r.CDRslow_mean <- l.inputs[["r.CDRslow_mean"]][t]
  } else {
    r.CDRslow_mean <- l.inputs[["r.CDRslow_mean"]]
  }
  
  # CDR (true)
  a.out[t,"CDR",alive] <- f.update_CDR(
    v.SYN            = a.out[t,"SYN",alive],
    v.SYN.lag        = a.out[t-1,"SYN",alive],
    v.MEMLOSS.lag    = a.out[t-1,"MEMLOSS",alive],
    cutoff_CDR       = l.inputs[["cutoff_CDR"]],
    v.CDR.lag        = a.out[t-1,"CDR",alive], 
    r.CDRfast_mean   = r.CDRfast_mean,
    r.CDRslow_mean   = r.CDRslow_mean,
    v.CDR_track      = a.out[t,"CDR_track",alive],
    v.CDRfast_sd1    = a.out[t,"CDRfast_sd1",alive],
    v.CDRslow_sd1    = a.out[t,"CDRslow_sd1",alive],
    r.CDR_sd2        = l.inputs[["r.CDR_sd2"]],
    v.TX.lag         = a.out[t-1,"TX",alive],
    rr.Tx_mci        = l.inputs[["rr.Tx_mci"]],
    random_cycle     = a.random[t,"CDR",alive],
    n.alive          = n.alive
  )
  
  # SEV (true)
  a.out[t,"SEV",alive] <- f.update_SEV(
    v.SYN          = a.out[t,"SYN",alive],
    v.CDR          = a.out[t,"CDR",alive],
    cutoff_CDR     = l.inputs[["cutoff_CDR"]],
    n.alive        = n.alive
  )
  
  return(a.out)
}



#### Module Functions ####
######################################## MEDBUR

f.update_MEDBUR <- function(v.MEDBUR.lag, v.AGE.lag, coef_MEDBUR, amplification, max_MEDBUR, random_cycle, n.alive) {
  
  # assume medbur distribution exp(0.08 * (age index)), scaled to 0-15
  medbur <- rep(NA, n.alive)
  
  prob <- 1 / (1 + exp(-coef_MEDBUR * (v.AGE.lag - 90)))
  prob_adj <- pmin((prob + amplification * v.MEDBUR.lag), 1)
  
  new_medbur <- round(qbinom(random_cycle, size = max_MEDBUR, prob = prob_adj), 0)
  
  medbur[new_medbur >= v.MEDBUR.lag] <- new_medbur[new_medbur >= v.MEDBUR.lag]
  medbur[new_medbur < v.MEDBUR.lag] <- v.MEDBUR.lag[new_medbur < v.MEDBUR.lag]
  
  return(medbur)
}

######################################## APOE4

f.update_APOE4 <- function(v.APOE4.lag) {
  apoe4 <- v.APOE4.lag
  return(apoe4)
}

######################################## SYN
f.update_SYN <- function(l.inputs, v.AGE.tplus2, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag, 
                         v.SYN.lag, v.SYN.lag2, v.MEMLOSS.lag,
                         random_tplus2, n.alive) {
  
  # start with empty vector
  symptoms <- rep(NA, n.alive)
  
  prob_mci <- f.calc_MCIprob(l.inputs, v.AGE.tplus2, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag)
  
  # apply dependent on previous state
  symptoms[v.SYN.lag == 0] <- as.numeric(prob_mci[v.SYN.lag == 0] > random_tplus2[v.SYN.lag == 0]) * 0.5
  symptoms[v.SYN.lag == 0.5 & v.SYN.lag2 == 0] <- 0.5
  symptoms[v.SYN.lag == 0.5 & v.SYN.lag2 == 0.5] <- 1
  symptoms[v.SYN.lag == 1 & v.MEMLOSS.lag == 1] <- as.numeric(0.07 < random_tplus2[v.SYN.lag == 1 & v.MEMLOSS.lag ==1])
  symptoms[v.SYN.lag == 1 & v.MEMLOSS.lag == 0] <- 1
  
  return(symptoms)
}

######################################## MEMLOSS - memory loss flag (reflects non-progressive impairment)

f.update_MEMLOSS <- function(v.MEMLOSS.lag, v.SYN, v.SYN.lag, v.AGE, v.EDU.lag, v.SEX.lag, v.RACEETH.lag, v.APOE4.lag, v.MEDBUR.lag, v.INCOME.lag, l.inputs, p.MEMLOSS_new, random_cycle, n.alive) {
  memloss <- rep(NA, n.alive)
  
  # New MCI cases: assign MEMLOSS with probability p.MEMLOSS_new
  new_mci <- v.SYN == 1 & v.SYN.lag < 1
  memloss[new_mci] <- as.numeric(p.MEMLOSS_new > random_cycle[new_mci])
  
  # For those with MEMLOSS==1, compute hazard for clearing (same as healthy-to-MCI hazard)
  idx_memloss <- which(v.MEMLOSS.lag == 1)
  if (length(idx_memloss) > 0) {
    prob_mci <- f.calc_MCIprob(l.inputs,
                               v.AGE[idx_memloss], v.EDU.lag[idx_memloss], v.SEX.lag[idx_memloss],
                               v.RACEETH.lag[idx_memloss], v.APOE4.lag[idx_memloss],
                               v.MEDBUR.lag[idx_memloss], v.INCOME.lag[idx_memloss])
    memloss_clear <- prob_mci > random_cycle[idx_memloss]
    memloss[idx_memloss[memloss_clear]] <- 0
    memloss[idx_memloss[!memloss_clear]] <- 1
  }
  
  # For those with MEMLOSS==0, keep at 0
  memloss[v.MEMLOSS.lag == 0 & !new_mci] <- 0
  
  # Ensure MEMLOSS is cleared if SYN == 0
  memloss[v.SYN == 0] <- 0
  
  return(memloss)
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

f.update_CDR <- function(v.SYN, v.SYN.lag, v.MEMLOSS.lag, cutoff_CDR, v.CDR.lag, 
                         r.CDRfast_mean, r.CDRslow_mean, v.CDR_track,
                         v.CDRfast_sd1, v.CDRslow_sd1, r.CDR_sd2, 
                         v.TX.lag, rr.Tx_mci, random_cycle, n.alive) {
  
  # CDR-SB cut-off values from: O'Bryant et al 2012 (PMC3409562) Table 2
  
  cdr <- rep(NA, n.alive)
  
  # Assign zero to those who are healthy
  healthy <- v.SYN == 0 & v.SYN.lag == 0
  cdr[healthy] <- 0
  
  re_healthy <- v.SYN == 0 & v.SYN.lag != 0
  cdr[re_healthy] <- qunif(p = random_cycle[re_healthy], min = cutoff_CDR["healthy"], max = cutoff_CDR["mci"])
  
  # Assign initial CDR-SB score for new cases (including MEMLOSS). Assume all enter in mild cognitive impairment
  mci_new <- v.SYN == 1 & v.SYN.lag < 1
  cdr[mci_new] <- qunif(p = random_cycle[mci_new], min = cutoff_CDR["mci"], max = cutoff_CDR["mild"])  # TODO: Skew right.
  
  # Keep CDR-SB score of those with prior memory loss (non-progressive impairment)
  mem_loss <- v.MEMLOSS.lag == 1 & !is.na(v.MEMLOSS.lag)
  cdr[mem_loss] <- v.CDR.lag[mem_loss]
  
  # Progress CDR-SB score for those already with impairment and no MEMLOSS
  mci_still <- v.SYN == 1 & v.SYN.lag == 1 & !mem_loss
  
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


