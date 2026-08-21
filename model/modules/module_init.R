######################################## GRAM MODULE: INITIALIZATION ########################################
# These functions initialize the model.

f.initialize <- function(l.inputs, microdata) {
  
  # generate enpty arrays for 1) random values and 2) outputs
  # with rows (1st dimension) = cycles, columns (2nd dimension) = attributes, 3rd dimension = individuals
  
  a.random <- a.out <- array(data = NA, dim = c(l.inputs[["n.cycle"]], l.inputs[["n.attr"]], l.inputs[["n.ind"]]),
                             dimnames = list(NULL, l.inputs[["v.attr_names"]], 1:l.inputs[["n.ind"]]))
  
  set.seed(l.inputs[["seed_stochastic"]])     # set the seed using the earlier defined seed
  a.random[,,] <- runif(n = length(a.random)) # put random value from uniform distribution in each element of the array
  
  n.ind <- l.inputs[["n.ind"]]
  
  # starting values (first cycle)
  a.out[1,"TIME",]      <- 0
  a.out[1,"ALIVE",]     <- 1
  
  if (!is.null(microdata)) {
    
    synthetic_pop <- generate_synthetic_sample(microdata, target_size = n.ind, weights = microdata$PERWT, seed = l.inputs[["seed_stochastic"]])
    
    a.out[1,"AGE",] <- synthetic_pop$AGE
    a.out[1,"SEX",] <- synthetic_pop$SEX
    a.out[1,"RACEETH",] <- synthetic_pop$RACEETH
    a.out[1,"EDU",] <- synthetic_pop$EDUC
    a.out[1,"INCOME",] <- synthetic_pop$INCOME_CAT
    a.out[1,"MEDBUR",] <- synthetic_pop$MEDBUR
    a.out[1,"APOE4",] <- synthetic_pop$APOE4
    a.out[1,"HCARE",] <- synthetic_pop$INSURANCE * as.numeric(a.random[1,"HCARE",] < l.inputs[["p.HCARE_start"]][2])
    
  } else {
    
    # Parametric cohort: every attribute is drawn from the distributions declared in setup.R,
    # so a cohort with no microdata (e.g. an external validation or country adaptation) is
    # specified entirely through l.inputs.
    
    # AGE is bounded by the age range the life table covers; ageing past that range is
    # handled by f.update_ALIVE, which holds mortality at the oldest tabulated age.
    AGE_min <- 50
    AGE_max <- AGE_min + nrow(l.inputs[["m.lifetable"]]) - 1
    a.out[1,"AGE",]       <- round(qnorm(p = a.random[1,"AGE",], mean = l.inputs[["AGE_start_mean"]], sd = l.inputs[["AGE_start_sd"]]),0)
    a.out[1,"AGE",]       <- pmin(pmax(a.out[1,"AGE",], AGE_min), AGE_max)
    
    a.out[1,"SEX",]       <- f.qcat(p_rand = a.random[1,"SEX",], p_cat = c(l.inputs[["p.SEX_start_male"]],
                                                                           l.inputs[["p.SEX_start_female"]]), values = l.inputs[["v.SEX_val"]])
    a.out[1,"RACEETH",]   <- f.qcat(p_rand = a.random[1,"RACEETH",], p_cat = l.inputs[["p.RACEETH_start"]], values = l.inputs[["v.RACEETH_val"]])
    
    # EDU (years). Drawn from the race-specific matrix m.EDU_start unless a cohort supplies a
    # single marginal distribution via p.EDU_start, or a continuous distribution via EDU_start_mean, which then take precedence.
    
    if (!is.null(l.inputs[["EDU_start_mean"]])) {
      a.out[1,"EDU",] <- round(qnorm(p = a.random[1,"EDU",], mean = l.inputs[["EDU_start_mean"]], sd = l.inputs[["EDU_start_sd"]]),0)
      a.out[1,"EDU",] <- pmin(pmax(a.out[1,"EDU",], 0), 26)
    } else if (!is.null(l.inputs[["p.EDU_start"]])) {
      a.out[1,"EDU",]     <- f.qcat(p_rand = a.random[1,"EDU",], p_cat = l.inputs[["p.EDU_start"]], values = l.inputs[["v.EDU_val"]])
    } else {
      v.EDU <- rep(NA, n.ind)
      for (i in seq_along(l.inputs[["v.RACEETH_val"]])) {
        idx <- which(a.out[1,"RACEETH",] == l.inputs[["v.RACEETH_val"]][i])
        if (length(idx) > 0) {
          v.EDU[idx] <- f.qcat(p_rand = a.random[1,"EDU",idx], p_cat = l.inputs[["m.EDU_start"]][ , i], values = l.inputs[["v.EDU_val"]])
        }
      }
      a.out[1,"EDU",]     <- v.EDU
    }
    
    a.out[1,"INCOME",]    <- f.qcat(p_rand = a.random[1,"INCOME",], p_cat = l.inputs[["p.INCOME_start"]], values = l.inputs[["v.INCOME_val"]])
    
    # MEDBUR: the chance of carrying 2+ conditions depends on education, and the number of
    # conditions is then drawn from the sex-specific prevalence of 0-10 conditions. Both stages
    # are folded into one per-individual probability matrix so a single random stream is used.
    
    if (!is.null(l.inputs[["MEDBUR_start_mean"]])) {
      a.out[1,"MEDBUR",] <- round(qnorm(p = a.random[1,"MEDBUR",], mean = l.inputs[["MEDBUR_start_mean"]], sd = l.inputs[["MEDBUR_start_sd"]]),0)
      a.out[1,"MEDBUR",] <- pmin(pmax(a.out[1,"MEDBUR",], 0), 15)
    } else {
      
      p.MEDBUR <- cbind(l.inputs[["p.MEDBUR_start_male"]], l.inputs[["p.MEDBUR_start_female"]])[ , a.out[1,"SEX",]]
      
      rr_EDU <- l.inputs[["rr.MEDBUR_2plus_EDU"]]
      p_2plus <- pmin(case_when(
        a.out[1,"EDU",] >= 16 ~ rr_EDU[["college"]],
        a.out[1,"EDU",] >= 12 ~ rr_EDU[["highschool"]],
        TRUE                  ~ rr_EDU[["lesshighschool"]]) * l.inputs[["p.MEDBUR_2plus_start"]], 1)
      
      m.MEDBUR <- rbind(
        sweep(p.MEDBUR[1:2, , drop = FALSE],  2, colSums(p.MEDBUR[1:2, , drop = FALSE]),  "/") * rep(1 - p_2plus, each = 2),
        sweep(p.MEDBUR[-(1:2), , drop = FALSE], 2, colSums(p.MEDBUR[-(1:2), , drop = FALSE]), "/") * rep(p_2plus, each = nrow(p.MEDBUR) - 2))
      
      a.out[1,"MEDBUR",]    <- f.qcat(p_rand = a.random[1,"MEDBUR",], p_cat = m.MEDBUR, values = seq_len(nrow(m.MEDBUR)) - 1)
    }
    
    a.out[1,"APOE4",]     <- f.qcat(p_rand = a.random[1,"APOE4",], p_cat = l.inputs[["p.APOE4_start"]], values = l.inputs[["v.APOE4_val"]])
    a.out[1,"HCARE",]     <- f.qcat(p_rand = a.random[1,"HCARE",], p_cat = l.inputs[["p.HCARE_start"]], values = l.inputs[["v.HCARE_val"]])
  }
  
  a.out[1,"TX",]         <- 0
  
  # SYN: 0 = healthy, 0.5 = TCI (a 2-cycle tunnel), 1 = impaired.
  a.out[1,"SYN",]        <- f.qcat(p_rand = a.random[1,"SYN",], p_cat = l.inputs[["p.SYN_start"]], values = l.inputs[["v.SYN_val"]])
  
  v.SYN_start <- a.out[1,"SYN",]
  tci_start   <- which(v.SYN_start == 0.5)
  impaired    <- which(v.SYN_start == 1)
  
  # TCI counts the cycles already spent in the TCI tunnel, and is 0 for anyone not in it.
  # Prevalent TCI cases have no history to read back, so their position in the tunnel is drawn.
  a.out[1,"TCI",]        <- 0
  if (length(tci_start) > 0) {
    a.out[1,"TCI",tci_start] <- f.qcat(p_rand = a.random[1,"TCI",tci_start], p_cat = l.inputs[["p.TCI_start"]], values = c(1,2))
  }
  
  a.out[1,"CDRfast_sd1",] <- qnorm(p = a.random[1,"CDRfast_sd1",], mean = 0, sd = l.inputs[["r.CDRfast_sd1"]])
  a.out[1,"CDRslow_sd1",] <- qnorm(p = a.random[1,"CDRslow_sd1",], mean = 0, sd = l.inputs[["r.CDRslow_sd1"]])
  
  # CDR-SB is 0 unless impaired. Prevalent impaired cases are assigned a severity from
  # p.SEV_start and then a CDR-SB score drawn uniformly within that severity's band, so that
  # CDR and SEV agree with each other and with the cutoffs used from cycle 2 onwards.
  cutoff <- l.inputs[["cutoff_CDR"]]
  a.out[1,"CDR",]        <- 0
  if (length(impaired) > 0) {
    v.SEV_start <- f.qcat(p_rand = a.random[1,"SEV",impaired], p_cat = l.inputs[["p.SEV_start"]], values = l.inputs[["v.SEV_val"]])
    CDR_lo <- cutoff[c("mci","mild","moderate","severe")][v.SEV_start + 1]
    CDR_hi <- cutoff[c("mild","moderate","severe","max")][v.SEV_start + 1]
    a.out[1,"CDR",impaired] <- qunif(p = a.random[1,"CDR",impaired], min = CDR_lo, max = CDR_hi)
  }
  
  # Derive SEV and CDR_track from CDR with the same functions used in every later cycle, so the
  # baseline obeys the same invariants (SEV is NA unless SYN == 1).
  a.out[1,"SEV",]        <- f.update_SEV(v.SYN = a.out[1,"SYN",], v.CDR = a.out[1,"CDR",],
                                         cutoff_CDR = cutoff, n.alive = n.ind)
  a.out[1,"CDR_track",]  <- f.update_CDR_track(v.SEV.lag = a.out[1,"SEV",], n.alive = n.ind)
  
  # MEMLOSS (non-progressive impairment) only applies to prevalent MCI cases.
  a.out[1,"MEMLOSS",]    <- 0
  mci_start <- which(a.out[1,"SYN",] == 1 & a.out[1,"SEV",] == 0)
  if (length(mci_start) > 0) {
    a.out[1,"MEMLOSS",mci_start] <- f.qcat(p_rand = a.random[1,"MEMLOSS",mci_start], p_cat = l.inputs[["p.MEMLOSS_start"]], values = l.inputs[["v.MEMLOSS_val"]])
  }
  
  a.out[1,"COGCON",]    <- 0    # assume no concerns at age 50 when healthy !! TODO: make this dynamic!
  a.out[1,"BHA",]       <- -9   # will be assigned as people are tested
  a.out[1,"last_BHA_age",]   <- NA   # tracking variable for implementing BHA scenarios
  a.out[1,"any_BHA_pos",]    <- FALSE
  a.out[1,"CDR_obs",]   <- -9   # DR-SB observed score -- will be assigned as people are tested
  a.out[1,"SEV_obs",]   <- -9   # will be assigned as people are tested
  
  a.out[1,"DX",]        <- f.qcat(p_rand = a.random[1,"DX",], p_cat = l.inputs[["p.DX_start"]], values = l.inputs[["v.DX_val"]])
  a.out[1,"PCP",]       <- NA
  a.out[1,"PET",]       <- NA
  a.out[1,"NP",]        <- NA
  a.out[1,"TX2",]       <- 0    # no one is on treatment at start
  a.out[1,"LTC",]       <- 0    # no one in long-term care at start
  a.out[1,"QALY",]      <- NA
  a.out[1,"COST_test",]  <- NA
  a.out[1,"COST_fu",]  <- NA
  a.out[1,"COST_tx2",]  <- NA
  a.out[1,"COST_care",] <- NA
  a.out[1,"COST_tx",]   <- NA
  
  return(list(
    a.random = a.random,
    a.out = a.out
  ))
}
