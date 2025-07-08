# ****************************************************************
# ======= 3. GRAM MAIN FUNCTIONS ======= 
# ****************************************************************

######################################## 3.1. MAIN FUNCTIONS: RUN MODEL ########################################

# run model
f.run <- function(l.inputs, microdata, printLevel) {
  
  # progress bar
  if (printLevel > 1) {
    ptm <- proc.time()
    stime <- Sys.time()
    pb = txtProgressBar(min = 0, max = l.inputs[["n.cycle"]], initial = 0)
  }
  
  # generate enpty arrays for 1) random values and 2) outputs
  # with rows (1st dimension) = cycles, columns (2nd dimension) = attributes, 3rd dimension = individuals
  a.random <- a.out <- array(data = NA, dim = c(l.inputs[["n.cycle"]], l.inputs[["n.attr"]], l.inputs[["n.ind"]]), 
                             dimnames = list(NULL, l.inputs[["v.attr_names"]], NULL))  
  
  set.seed(l.inputs[["seed_stochastic"]])     # set the seed using the earlier defined seed
  a.random[,,] <- runif(n = length(a.random)) # put random value from uniform distribution in each element of the array
  
  # starting values (first cycle)
  a.out[1,"TIME",]      <- 0
  a.out[1,"ALIVE",]     <- 1
  
  if (!is.null(microdata)) {

    synthetic_pop <- generate_synthetic_sample(microdata, target_size = l.inputs[["n.ind"]], seed = l.inputs[["seed_stochastic"]])

    a.out[1,"AGE",] <- synthetic_pop$AGE
    a.out[1,"SEX",] <- synthetic_pop$SEX
    a.out[1,"RACEETH",] <- synthetic_pop$RACEETH
    a.out[1,"EDU",] <- synthetic_pop$EDUC
    a.out[1,"INCOME",] <- synthetic_pop$INCOME_CAT
    a.out[1,"MEDBUR",] <- synthetic_pop$MEDBUR
    a.out[1,"APOE4",] <- synthetic_pop$APOE4
    a.out[1,"HCARE",] <- synthetic_pop$INSURANCE * as.numeric(a.random[1,"HCARE",] < l.inputs[["p.HCARE_start"]][2]) 
  } else {
    a.out[1,"AGE",]       <- round(qnorm(p = a.random[1,"AGE",], mean = l.inputs[["AGE_start_mean"]], sd = l.inputs[["AGE_start_sd"]]),0)
    a.out[1,"AGE",][a.out[1,"AGE",]<50] <- 50
    a.out[1,"AGE",][a.out[1,"AGE",]>99] <- 99
    a.out[1,"SEX",]       <- f.qcat(p_rand = a.random[1,"SEX",], p_cat = c(l.inputs[["p.SEX_start_male"]], 
                                                                           l.inputs[["p.SEX_start_female"]]), values = l.inputs[["v.SEX_val"]])
    a.out[1,"RACEETH",]   <- f.qcat(p_rand = a.random[1,"RACEETH",], p_cat = l.inputs[["p.RACEETH_start"]], values = l.inputs[["v.RACEETH_val"]])
    
    a.out[1,"EDU",]       <- case_match(a.out[1,"RACEETH",],
                                        0 ~ f.qcat(p_rand = a.random[1,"EDU",], 
                                                   p_cat = c(0.048, 0.274, 0.149, 0.111, 0.261, 0.117, 0.017, 0.023),
                                                   values = c(8, 12, 14, 14, 16, 18, 19, 25)),
                                        1 ~ f.qcat(p_rand = a.random[1,"EDU",], 
                                                   p_cat = c(0.095, 0.335, 0.181, 0.110, 0.173, 0.081, 0.010, 0.015),
                                                   values = c(8, 12, 14, 14, 16, 18, 19, 25)),
                                        2 ~ f.qcat(p_rand = a.random[1,"EDU",], 
                                                   p_cat = c(0.248, 0.327, 0.130, 0.086, 0.145, 0.047, 0.009, 0.008),
                                                   values = c(8, 12, 14, 14, 16, 18, 19, 25)))
    # probs from https://www.equityinhighered.org/indicators/u-s-population-trends-and-educational-attainment/educational-attainment-by-race-and-ethnicity/
    #  round(qbeta(p = a.random[1,"EDU",], shape1 = 2.24, shape2 = 2.90) * 25, 0)  # assumes max 25 years of education
    a.out[1,"INCOME",]    <- f.qcat(p_rand = a.random[1,"INCOME",], p_cat = l.inputs[["p.INCOME_start"]], values = l.inputs[["v.INCOME_val"]])
    
    tmp_rand <- runif(n = l.inputs[["n.ind"]])
    tmp_probs_2plus_medbur <- case_when(
      a.out[1,"EDU",] >= 16 ~ 1,
      a.out[1,"EDU",] >= 12 ~ 1.32,
      a.out[1,"EDU",] < 12 ~ 1.58) * 0.531    # See script "calculate_initial_medbur.R"
    tmp_2plus_medbur <- qbinom(p = tmp_rand, size = 1, prob = tmp_probs_2plus_medbur)
    
    prev.no_of_conditions_male <- c(0.155, 0.205, 0.220, 0.175, 0.105, 0.060, 0.045, 0.020, 0.018, 0.005, 0.001)
    prev.no_of_conditions_female <- c(0.175, 0.205, 0.190, 0.175, 0.140, 0.075, 0.040, 0.018, 0.015, 0.002, 0.002)
    prev.no_of_conditions <- (prev.no_of_conditions_female + prev.no_of_conditions_male) / 2
    
    a.out[1,"MEDBUR",]    <- case_when(
      tmp_2plus_medbur == 1 ~ f.qcat(p_rand = a.random[1,"MEDBUR",], p_cat = prev.no_of_conditions[3:11]/sum(prev.no_of_conditions[3:11]), values = 2:10),
      tmp_2plus_medbur == 0 ~ f.qcat(p_rand = a.random[1,"MEDBUR",], p_cat = prev.no_of_conditions[1:2]/sum(prev.no_of_conditions[1:2]), values = 0:1)
    )
    
    # round(qbeta(p = a.random[1,"MEDBUR",], shape1 = 2, shape2 = 18) * max(l.inputs[["v.MEDBUR_val"]]),0)
    a.out[1,"APOE4",]     <- f.qcat(p_rand = a.random[1,"APOE4",], p_cat = l.inputs[["p.APOE4_start"]], values = l.inputs[["v.APOE4_val"]])
    a.out[1,"HCARE",]     <- f.qcat(p_rand = a.random[1,"HCARE",], p_cat = l.inputs[["p.HCARE_start"]], values = l.inputs[["v.HCARE_val"]])
  }
  
  a.out[1,"TX",]         <- 0
  
  a.out[1,"SYN",]        <- 0   # everyone starts healthy
  
  a.out[1,"CDR_track",]  <- case_match(a.out[1,"SYN",],
                                       0 ~ 0,
                                       1 ~ 1)
  a.out[1,"CDRfast_sd1",] <- qnorm(p = a.random[1,"CDRfast_sd1",], mean = 0, sd = l.inputs[["r.CDRfast_sd1"]]) 
  a.out[1,"CDRslow_sd1",] <- qnorm(p = a.random[1,"CDRslow_sd1",], mean = 0, sd = l.inputs[["r.CDRslow_sd1"]]) 
  a.out[1,"CDR",]       <- case_match(a.out[1,"SYN",],
                                      0 ~ 0,
                                      1 ~ qunif(p = a.random[1,"CDR",], min = l.inputs[["cutoff_CDR"]]["mci"], max = l.inputs[["cutoff_CDR"]]["moderate"]))
  a.out[1,"MEMLOSS",]   <- case_match(a.out[1,"SYN",],
                                      0 ~ 0,
                                      1 ~ f.qcat(p_rand = a.random[1,"MEMLOSS",], p_cat = l.inputs[["p.MEMLOSS_start"]], values = l.inputs[["v.MEMLOSS_val"]]))   # !! when starting with prevalent pop, this will cause a problem as people with dem might get MEMLOSS of 1
  a.out[1,"SEV",]       <- f.update_SEV(v.SYN = a.out[1,"SYN",], v.CDR = a.out[1,"CDR",], 
                                        cutoff_CDR = l.inputs[["cutoff_CDR"]], n.alive = l.inputs[["n.ind"]])   # Everyone starts healthy
  a.out[1,"COGCON",]    <- 0    # Assume no concerns at age 50 when healthy !! TODO: make this dynamic!
  a.out[1,"BHA",]       <- -9   # Will be assigned as people are tested
  a.out[1,"CDR_obs",]   <- -9   # CDR-SB observed score -- will be assigned as people are tested
  a.out[1,"SEV_obs",]   <- -9   # Will be assigned as people are tested
  
  a.out[1,"DX",]        <- f.qcat(p_rand = a.random[1,"DX",], p_cat = l.inputs[["p.DX_start"]], values = l.inputs[["v.DX_val"]])
  a.out[1,"FUN",]       <- 0    # FAQ -- everyone starts with no functional impairment
  a.out[1,"BEH",]       <- NA   # Currently not in use
  a.out[1,"INSTIT",]    <- 0    # Currently not in use
  a.out[1,"QALY",]      <- NA
  a.out[1,"COST_care",] <- NA
  a.out[1,"COST_tx",]   <- NA
  
  # update progress bar
  if (printLevel > 1) setTxtProgressBar(pb, 1)
  
  # run subsequent cycles
  for(t in 2:l.inputs[["n.cycle"]]) {
    
    # TIME
    a.out[t,"TIME",] <- f.update_TIME(
      v.TIME.lag   = a.out[t-1,"TIME",]
    )
    
    alive.lag <- a.out[t-1,"ALIVE",]==1
    
    # ALIVE
    a.out[t,"ALIVE",alive.lag] <- f.update_ALIVE(
      alive.lag     = alive.lag,
      v.AGE.lag     = a.out[t-1,"AGE",alive.lag], 
      v.SYN.lag     = a.out[t-1,"SYN",alive.lag],
      v.SEV.lag     = a.out[t-1,"SEV",alive.lag], 
      random_cycle  = a.random[t,"ALIVE",alive.lag], 
      m.lifetable   = l.inputs[["m.lifetable"]], 
      hr.mort_mci   = l.inputs[["hr.mort_mci"]], 
      hr.mort_mil   = l.inputs[["hr.mort_mil"]], 
      hr.mort_mod   = l.inputs[["hr.mort_mod"]], 
      hr.mort_sev   = l.inputs[["hr.mort_sev"]],
      hr.mort_mci_age   = l.inputs[["hr.mort_mci_age"]],
      hr.mort_mil_age   = l.inputs[["hr.mort_mil_age"]],
      hr.mort_mod_age   = l.inputs[["hr.mort_mod_age"]],
      hr.mort_sev_age   = l.inputs[["hr.mort_sev_age"]]
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
    
    # HCARE
    a.out[t,"HCARE",alive] <- f.update_HCARE(
      v.HCARE.lag  = a.out[t-1,"HCARE",alive],
      v.AGE        = a.out[t,"AGE",alive],
      random_cycle = a.random[t,"HCARE",alive]
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
    
    # COGCON
    a.out[t,"COGCON",alive] <- f.update_COGCON(
      v.AGE            = a.out[t,"AGE",alive],
      v.SYN            = a.out[t,"SYN",alive],
      v.SEV            = a.out[t,"SEV",alive],
      m.cogcon         = l.inputs[["m.cogcon"]],
      v.DX.lag         = a.out[t-1,"DX",alive],
      random_cycle     = a.random[t,"COGCON",alive],
      n.alive          = n.alive
    )
    
    # BHA
    a.out[t,"BHA",alive] <- f.update_BHA(
      v.HCARE          = a.out[t,"HCARE",alive],
      v.DX.lag         = a.out[t-1,"DX",alive],
      v.COGCON         = a.out[t,"COGCON",alive],
      v.SYN            = a.out[t,"SYN",alive],
      v.SEV            = a.out[t,"SEV",alive],
      v.MEMLOSS        = a.out[t,"MEMLOSS",alive], 
      sens_BHA         = l.inputs[["sens_BHA"]],
      spec_BHA         = l.inputs[["spec_BHA"]],
      random_cycle     = a.random[t,"BHA",alive],
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
    
    # SEV_obs
    a.out[t,"SEV_obs",alive] <- f.update_SEV_obs(
      v.CDR_obs     = a.out[t,"CDR_obs",alive],
      cutoff_CDR     = l.inputs[["cutoff_CDR"]],
      n.alive        = n.alive
    ) 
    
    
    # DX
    a.out[t,"DX",alive] <- f.update_DX(
      v.DX.lag       = a.out[t-1,"DX",alive],
      v.SYN          = a.out[t,"SYN",alive], 
      v.SEV          = a.out[t,"SEV",alive], 
      v.HCARE        = a.out[t,"HCARE",alive],
      random_cycle   = a.random[t,"DX",alive], 
      n.alive        = n.alive
    ) 
    
    
    # FAQ 
    
    
    
    # TX
    
    Tx_counter <- a.out[t-1,"TX",]
    
    a.out[t,"TX",alive] <- f.update_TX(
      v.SEV_obs            = a.out[t,"SEV_obs",alive], 
      v.SEV_obs.lag        = a.out[t-1,"SEV_obs",alive],
      Tx_counter           = Tx_counter[alive],
      v.TX.lag             = a.out[t-1,"TX",alive], 
      n.alive              = n.alive, 
      random_cycle         = a.random[t,"TX",alive], 
      Tx_t_max             = l.inputs[["Tx_t_max"]], 
      p.Tx                 = l.inputs[["p.Tx"]], 
      Tx                   = l.inputs[["Tx"]]
    )
    
    Tx_counter[alive] <- Tx_counter[alive] + a.out[t,"TX",alive]
    
    
    # update progress bar
    if (printLevel > 1) setTxtProgressBar(pb, t)
    
  }
  
  # Convert all non-observed values to NA
  # a.out[a.out == -9] <- NA
  
  # The "round(score*2)/2" ensures that scores are increments of 0.5 
  a.out[,"CDR",] <- round(a.out[,"CDR",] * 2) / 2
  a.out[,"CDR_obs",] <- round(a.out[,"CDR_obs",] * 2) / 2
  
  
  
  # run time
  if (printLevel >1 ) {
    print(proc.time() - ptm)
    print(Sys.time() - stime)
  }
  
  # return outcome
  return(a.out)
  
}

######################################## 3.2. MAIN FUNCTIONS: GENERATE RESULTS ########################################


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
  
  n <- dim(a.out)[3]     # get the number of individuals in the subset
  
  # mean time alive
  l.out[["alive.trc"]] <- as.matrix(apply(X = a.out[,"ALIVE",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["alive.sum"]] <- sum(l.out[["alive.trc"]])
  l.out[["alive.trc.dis"]] <- as.matrix(f.discount(x = l.out[["alive.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["alive.sum.dis"]] <- sum(l.out[["alive.trc.dis"]])
  
  # mean time healthy
  l.out[["healthy.trc"]] <- as.matrix(apply(X = a.out[,"SYN",]<1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["healthy.sum"]] <- sum(l.out[["healthy.trc"]])
  l.out[["healthy.trc.dis"]] <- as.matrix(f.discount(x = l.out[["healthy.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["healthy.sum.dis"]] <- sum(l.out[["healthy.trc.dis"]])
  l.out[["healthy_obs.trc"]] <- as.matrix(apply(X = a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["healthy_obs.sum"]] <- sum(l.out[["healthy_obs.trc"]])
  
  # time in MCI
  l.out[["MCI.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["MCI.sum"]] <- sum(l.out[["MCI.trc"]])
  l.out[["MCI.trc.dis"]] <- as.matrix(f.discount(x = l.out[["MCI.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["MCI.sum.dis"]] <- sum(l.out[["MCI.trc.dis"]])
  l.out[["MCI_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["MCI_obs.sum"]] <- sum(l.out[["MCI_obs.trc"]])
  
  # time in mild dementia
  l.out[["SEV1.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV1.sum"]] <- sum(l.out[["SEV1.trc"]])
  l.out[["SEV1.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV1.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV1.sum.dis"]] <- sum(l.out[["SEV1.trc.dis"]])
  l.out[["SEV1_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV1_obs.sum"]] <- sum(l.out[["SEV1_obs.trc"]])
  
  # time in moderate dementia
  l.out[["SEV2.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV2.sum"]] <- sum(l.out[["SEV2.trc"]])
  l.out[["SEV2.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV2.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV2.sum.dis"]] <- sum(l.out[["SEV2.trc.dis"]])
  l.out[["SEV2_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV2_obs.sum"]] <- sum(l.out[["SEV2_obs.trc"]])
  
  # time in severe dementia
  l.out[["SEV3.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV3.sum"]] <- sum(l.out[["SEV3.trc"]])
  l.out[["SEV3.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV3.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV3.sum.dis"]] <- sum(l.out[["SEV3.trc.dis"]])
  l.out[["SEV3_obs.trc"]] <- as.matrix(apply(X = a.out[,"SEV_obs",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["SEV3_obs.sum"]] <- sum(l.out[["SEV3_obs.trc"]])
  
  # time in treatment
  l.out[["TX.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["TX.sum"]] <- sum(sum(l.out[["TX.trc"]]))
  
  # temporary to select outcomes
  l.out[["mean_time_alive"]] <- l.out[["alive.sum.dis"]]
  l.out[["mean_time_healthy"]] <- l.out[["healthy.sum.dis"]]
  l.out[["mean_time_MCI"]] <- l.out[["MCI.sum.dis"]]
  l.out[["mean_time_SEV1"]] <- l.out[["SEV1.sum.dis"]]
  l.out[["mean_time_SEV2"]] <- l.out[["SEV2.sum.dis"]]
  l.out[["mean_time_SEV3"]] <- l.out[["SEV3.sum.dis"]]
  
  # time in treatment
  l.out[["mean_time_Tx"]] <- sum(a.out[,"TX",]==1, na.rm=TRUE)/n
  
  # age at onset
  onset_cycle <- apply(a.out[,"SYN",], 2, function(x) {
    # Identify indices where SYN == 1
    mci_indices <- which(!is.na(x) & x == 1)
    
    # Return the first such index, or NA if none exist
    if (length(mci_indices) > 0) mci_indices[1] else NA
  })
  
  l.out[["age_at_onset"]] <- ifelse(!is.na(onset_cycle), a.out[cbind(onset_cycle, match("AGE", dimnames(a.out)[[2]]), seq_len(dim(a.out)[3]))], NA)
  
  # state reside time by age at onset
  
  last_alive_cycle <- apply(a.out[, "ALIVE", ], 2, function(x) which(x == 0)[1])-1
  sev_at_death <- ifelse(!is.na(last_alive_cycle), a.out[cbind(last_alive_cycle, match("SEV", dimnames(a.out)[[2]]), seq_len(dim(a.out)[3]))], NA)
  
  age_bins <- cut(l.out[["age_at_onset"]], breaks = seq(50, 100, by = 5), right = FALSE, include.lowest = TRUE)
  age_groups <- levels(age_bins)
  severity_levels <- c("mci", "mil", "mod", "sev")
  result_matrix1 <- result_matrix2 <- matrix(NA, nrow = length(age_groups)+1, ncol = length(severity_levels),
                                             dimnames = list(c(age_groups,"Overall"), severity_levels))
  
  for (sev in 0:3) {
    # Total time in the current state for each individual
    time_in_state <- colSums(a.out[,"SEV",] == sev, na.rm = TRUE)
    
    if (sev < 3) {
      time_in_state_censored <- ifelse(sev_at_death == sev, NA, time_in_state)
    } else {
      time_in_state_censored <- time_in_state
    }
    
    # Average time per age group
    result_matrix1[1:length(age_groups), sev + 1] <- tapply(time_in_state, age_bins, function(x) round(mean(x[x > 0], na.rm = TRUE), digits = 2))
    result_matrix2[1:length(age_groups), sev + 1] <- tapply(time_in_state_censored, age_bins, function(x) round(mean(x[x > 0], na.rm = TRUE), digits = 2))
    
    # Average time overall
    result_matrix1[nrow(result_matrix1), sev + 1] <- round(mean(time_in_state[time_in_state > 0], na.rm = TRUE), digits = 2)
    result_matrix2[nrow(result_matrix2), sev + 1] <- round(mean(time_in_state_censored[time_in_state_censored > 0], na.rm = TRUE), digits = 2)

  }
  
  time_in_dem <- colSums(a.out[,"SEV",] > 0, na.rm = TRUE)
  time_in_dem_grouped <- tapply(time_in_dem, age_bins, function(x) round(mean(x[x > 0], na.rm = TRUE), digits = 2))
  time_in_dem_overall <- round(mean(time_in_dem[time_in_dem > 0], na.rm = TRUE), digits = 2)
  
  
  result_df1 <- as.data.frame(result_matrix1) %>% 
    mutate(age_group = factor(c(age_groups, "Overall"),
                              labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+", "Overall"))) %>%
    select(age_group, mci, mil, mod, sev) %>%
    mutate(any_dem = c(time_in_dem_grouped, time_in_dem_overall))
  rownames(result_df1) <- NULL
  
  result_df2 <- as.data.frame(result_matrix2) %>% 
    mutate(age_group = factor(c(age_groups, "Overall"),
                              labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+", "Overall"))) %>%
    select(age_group, mci, mil, mod, sev) %>%
    mutate(any_dem = c(time_in_dem_grouped, time_in_dem_overall))
  rownames(result_df2) <- NULL
  
  l.out[["reside_time"]] <- list(noncensored = result_df1, censored = result_df2)
  
  # time in full-time care
  # l.out[["mean_time_FTC"]] <- sum(a.out[,"INSTIT",]==1, na.rm=TRUE)/n
  
  # state trace (undiscounted) (true states)
  l.out[["state_trace"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 6, 
                                   dimnames = list(NULL,c("healthy","mci","mil","mod","sev","dth")))
  l.out[["state_trace"]][,"healthy"] <- l.out[["healthy.trc"]]
  l.out[["state_trace"]][,"mci"] <- l.out[["MCI.trc"]]
  l.out[["state_trace"]][,"mil"] <- l.out[["SEV1.trc"]]
  l.out[["state_trace"]][,"mod"] <- l.out[["SEV2.trc"]]
  l.out[["state_trace"]][,"sev"] <- l.out[["SEV3.trc"]]
  l.out[["state_trace"]][,"dth"] <- apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n
  # check rowsum
  rowSums(l.out[["state_trace"]])
  # trace institutionalized
  # l.out[["state_trace_instit"]] <- as.matrix((apply(X = a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n))
  
  # state trace (undiscounted) (observed states)
  l.out[["state_trace_obs"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 6, 
                                       dimnames = list(NULL,c("healthy","mci","mil","mod","sev","dth")))
  l.out[["state_trace_obs"]][,"healthy"] <- l.out[["healthy_obs.trc"]]
  l.out[["state_trace_obs"]][,"mci"] <- l.out[["MCI_obs.trc"]]
  l.out[["state_trace_obs"]][,"mil"] <- l.out[["SEV1_obs.trc"]]
  l.out[["state_trace_obs"]][,"mod"] <- l.out[["SEV2_obs.trc"]]
  l.out[["state_trace_obs"]][,"sev"] <- l.out[["SEV3_obs.trc"]]
  l.out[["state_trace_obs"]][,"dth"] <- apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n
  # check rowsum
  rowSums(l.out[["state_trace_obs"]])
  
  # correct/incorrect detection
  l.out[["state_concordance"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 13,
                                         dimnames = list(NULL, c("h_NA","h_neg","h_pos",
                                                                 "tci_NA", "tci_neg","tci_pos",
                                                                 "mci_NA","mci_neg","mci_pos",
                                                                 "dem_NA","dem_neg","dem_pos",
                                                                 "dth")))
  l.out[["state_concordance"]][,"h_NA"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==-9, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"h_neg"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"h_pos"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"tci_NA"] <- as.matrix(apply(X = a.out[,"SYN",]==0.5 & a.out[,"BHA",]==-9, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"tci_neg"] <- as.matrix(apply(X = a.out[,"SYN",]==0.5 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"tci_pos"] <- as.matrix(apply(X = a.out[,"SYN",]==0.5 & a.out[,"BHA",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_NA"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==-9, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_neg"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_pos"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_NA"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==-9, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_neg"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_pos"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dth"] <- as.matrix(apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  rowSums(l.out[["state_concordance"]])
  
  
  # MCI incidence
  temp.age_groups <- seq(50, 100, by= 5)
  l.out[["mci_incidence"]] <- matrix(NA, nrow = length(temp.age_groups), ncol = 1, 
                                     dimnames = list(paste0(temp.age_groups, "-", temp.age_groups + 4), "Incidence"))
  
  for (i in seq_along(temp.age_groups)) {
    age_min <- temp.age_groups[i]
    age_max <- age_min + 4
    new_cases <- 0
    person_years <- 0
    
    for (t in 2:dim(a.out)[1]) {  # Start from 2 to access previous cycle
      # Select individuals in the age group who were healthy at the start of the cycle
      at_risk <- which(a.out[t - 1, "AGE", ] >= age_min & a.out[t - 1, "AGE", ] <= age_max & a.out[t - 1, "SYN", ] != 1)
      
      # Calculate person-years at risk for this time point
      person_years <- person_years + length(at_risk)
      
      # Count new MCI cases (healthy in previous cycle and MCI in current cycle)
      new_cases <- new_cases + sum(a.out[t, "SYN", at_risk] == 1, na.rm = TRUE)
    }
    
    incidence_rate <- if (person_years > 0) new_cases / person_years else NA
    l.out[["mci_incidence"]][i, 1] <- incidence_rate
  }
  
  # Prevalence by severity by age -- this requires that age and cycle are equal (everyone starts the same age)
  l.out[["prevalence_by_age"]] <- as.data.frame(l.out[["state_trace"]]) %>%
    mutate(age = l.inputs[["AGE_start_mean"]]:(l.inputs[["AGE_start_mean"]]+l.inputs[["n.cycle"]]-1)) %>%
    mutate(alive = 1 - dth,
           age_group = factor(cut(age, breaks = c(50, 65, 75, 85, Inf), right = FALSE),
                              labels = c("50-64", "65-74","75-84","85+")
           )) %>%
    mutate(mci = mci / alive,
           mil = mil / alive,
           mod = mod / alive,
           sev = sev / alive) %>%
    group_by(age_group) %>%
    summarize(avg_mci = mean(mci),
              avg_mil = mean(mil),
              avg_mod = mean(mod),
              avg_sev = mean(sev)) %>%
    bind_rows(
      summarize(., 
                age_group = "Overall",
                avg_mci = mean(avg_mci), 
                avg_mil = mean(avg_mil), 
                avg_mod = mean(avg_mod), 
                avg_sev = mean(avg_sev))
    )
  
  # Prevalence by race/ethnicity
  l.out[["prevalence_by_raceeth"]] <- array(NA, dim = c(l.inputs[["n.cycle"]], length(severity_levels) + 1, length(l.inputs[["v.RACEETH_val"]])),
                                            dimnames = list(NULL, c("h", severity_levels), c("NHW","NHB","Hisp")))
  
  # matrix_nhw <- matrix_nhb <- matrix_hisp <- matrix(0, nrow = l.inputs[["n.cycle"]], ncol = 5,
  #                     dimnames = list(NULL, c("h", severity_levels)))
  l.out[["prevalence_by_raceeth"]][,"h","NHW"]   <- as.matrix(apply(X = a.out[,"RACEETH",]==0 & a.out[,"SYN",]<1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mci","NHW"] <- as.matrix(apply(X = a.out[,"RACEETH",]==0 & a.out[,"SEV",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mil","NHW"] <- as.matrix(apply(X = a.out[,"RACEETH",]==0 & a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mod","NHW"] <- as.matrix(apply(X = a.out[,"RACEETH",]==0 & a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"sev","NHW"] <- as.matrix(apply(X = a.out[,"RACEETH",]==0 & a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,,"NHW"] <- round(l.out[["prevalence_by_raceeth"]][,,"NHW"] / rowSums(l.out[["prevalence_by_raceeth"]][,,"NHW"]), 3)
  
  
  l.out[["prevalence_by_raceeth"]][,"h","NHB"]   <- as.matrix(apply(X = a.out[,"RACEETH",]==1 & a.out[,"SYN",]<1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mci","NHB"] <- as.matrix(apply(X = a.out[,"RACEETH",]==1 & a.out[,"SEV",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mil","NHB"] <- as.matrix(apply(X = a.out[,"RACEETH",]==1 & a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mod","NHB"] <- as.matrix(apply(X = a.out[,"RACEETH",]==1 & a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"sev","NHB"] <- as.matrix(apply(X = a.out[,"RACEETH",]==1 & a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,,"NHB"] <- round(l.out[["prevalence_by_raceeth"]][,,"NHB"] / rowSums(l.out[["prevalence_by_raceeth"]][,,"NHB"]), 3)
  
  l.out[["prevalence_by_raceeth"]][,"h","Hisp"]   <- as.matrix(apply(X = a.out[,"RACEETH",]==2 & a.out[,"SYN",]<1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mci","Hisp"] <- as.matrix(apply(X = a.out[,"RACEETH",]==2 & a.out[,"SEV",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mil","Hisp"] <- as.matrix(apply(X = a.out[,"RACEETH",]==2 & a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"mod","Hisp"] <- as.matrix(apply(X = a.out[,"RACEETH",]==2 & a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,"sev","Hisp"] <- as.matrix(apply(X = a.out[,"RACEETH",]==2 & a.out[,"SEV",]==3, MARGIN = 1, FUN = sum, na.rm = TRUE))
  l.out[["prevalence_by_raceeth"]][,,"Hisp"] <- round(l.out[["prevalence_by_raceeth"]][,,"Hisp"] / rowSums(l.out[["prevalence_by_raceeth"]][,,"Hisp"]), 3)
  
  
  
  # QALY
  l.out[["QALY"]] <- as.matrix(apply(X = a.out[,"QALY",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["QALY.sum"]] <- sum(l.out[["QALY"]])
  l.out[["QALY.dis"]] <- as.matrix(f.discount(x = l.out[["QALY"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["QALY.dis.sum"]] <- sum(l.out[["QALY.dis"]])
  
  # COST_test
  # l.out[["COST_test"]] <- as.matrix(apply(X = a.out[,"COST_test",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  # l.out[["COST_test.sum"]] <- sum(l.out[["COST_test"]])
  # l.out[["COST_test.dis"]] <- as.matrix(f.discount(x = l.out[["COST_test"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  # l.out[["COST_test.dis.sum"]] <- sum(l.out[["COST_test.dis"]])
  
  # COST_tx
  l.out[["COST_tx"]] <- as.matrix(apply(X = a.out[,"COST_tx",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["COST_tx.sum"]] <- sum(l.out[["COST_tx"]])
  l.out[["COST_tx.dis"]] <- as.matrix(f.discount(x = l.out[["COST_tx"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_tx.dis.sum"]] <- sum(l.out[["COST_tx.dis"]])
  
  # COST_fu
  # l.out[["COST_fu"]] <- as.matrix(apply(X = a.out[,"COST_fu",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  # l.out[["COST_fu.sum"]] <- sum(l.out[["COST_fu"]])
  # l.out[["COST_fu.dis"]] <- as.matrix(f.discount(x = l.out[["COST_fu"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  # l.out[["COST_fu.dis.sum"]] <- sum(l.out[["COST_fu.dis"]])
  # 
  # COST_care
  l.out[["COST_care"]] <- as.matrix(apply(X = a.out[,"COST_care",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
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
f.out_summary <- function(l.inputs1, l.inputs2, printLevel = 0) {
  
  # strategy 1
  l.inputs_strat1 <- l.inputs1
  # l.inputs_strat1[["strategy"]] <- l.inputs[["strategy_strat1"]]
  # l.inputs_strat1[["Tx"]] <- l.inputs[["Tx_strat1"]]
  a.out_strat1 <- f.run(l.inputs = l.inputs_strat1, printLevel = printLevel)
  a.out_qc_strat1 <- f.qaly_cost(a.out = a.out_strat1, l.inputs = l.inputs_strat1)
  out_strat1 <- f.out_aggregate(a.out = a.out_qc_strat1, l.inputs = l.inputs_strat1)
  fig_strat1 <- f.make_figures(l.out = out_strat1, l.inputs = l.inputs1)
  
  
  # strategy 2
  l.inputs_strat2 <- l.inputs2
  # l.inputs_strat2[["strategy"]] <- l.inputs[["strategy_strat2"]]
  # l.inputs_strat2[["Tx"]] <- l.inputs[["Tx_strat2"]]
  a.out_strat2 <- f.run(l.inputs = l.inputs_strat2, printLevel = printLevel)
  a.out_qc_strat2 <- f.qaly_cost(a.out = a.out_strat2, l.inputs = l.inputs_strat2)
  out_strat2 <- f.out_aggregate(a.out = a.out_qc_strat2, l.inputs = l.inputs_strat2)
  fig_strat2 <- f.make_figures(l.out = out_strat2, l.inputs = l.inputs2)
  
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
    l.inputs1 = l.inputs1,
    l.inputs2 = l.inputs2,
    a.out_strat1 = a.out_strat1, a.out_strat2 = a.out_strat2,
    a.out_qc_strat1 = a.out_qc_strat1, a.out_qc_strat2 = a.out_qc_strat2,
    out_strat1 = out_strat1, out_strat2 = out_strat2,
    fig_strat1 = fig_strat1, fig_strat2 = fig_strat2,
    m.out = m.out
  ))
  
}



######################################## 3.3. MAIN FUNCTIONS: WRAPPER  ########################################

f.wrap_run <- function(l.inputs, microdata = NULL, printLevel = 0) {
  
  
  output <- f.run(l.inputs = l.inputs, microdata = microdata, printLevel = printLevel)
  
  # Aggregate results for full cohort
  aggregated_results_totpop <- f.out_aggregate(a.out = output, l.inputs = l.inputs)
  fig.progression <- f.make_figures(l.out = aggregated_results_totpop, l.inputs = l.inputs)
  
  # # Aggregate results for those who develop MCI at any point
  # impaired <- apply(output[,"SYN",], 2, function(x) any(x == 1, na.rm = TRUE))
  # output_impaired <- output[,,impaired]
  # aggregated_results_impaired <- f.out_aggregate(a.out = output_impaired, l.inputs = l.inputs)
  # figures_impaired <- f.figures(l.out = aggregated_results_impaired, l.inputs = l.inputs)
  # 
  # # Aggregate results for only DMT eligibles
  # treated <- apply(output[,"TX",], 2, function(x) any(x == 1, na.rm = TRUE))
  # output_treated <- output[,,treated]
  # aggregated_results_treated <- f.out_aggregate(a.out = output_treated, l.inputs = l.inputs)
  # figures_treated <- f.figures(l.out = aggregated_results_treated, l.inputs = l.inputs)
  # 
  return(list(
    output = output,
    aggregated_results_totpop = aggregated_results_totpop,
    fig.progression = fig.progression
    # ,
    # aggregated_results_impaired = aggregated_results_impaired,
    # figures_impaired = figures_impaired,
    # aggregated_results_treated = aggregated_results_treated,
    # figures_treated = figures_treated
  )
  )
}

