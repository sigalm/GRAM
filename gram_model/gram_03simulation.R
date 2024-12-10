# ****************************************************************
# ======= 3. GRAM MAIN FUNCTIONS ======= 
# ****************************************************************

######################################## 3.1. MAIN FUNCTIONS: RUN MODEL ########################################

# run model
f.run <- function(l.inputs, printLevel) {
  
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
  a.out[1,"AGE",]       <- round(qnorm(p = a.random[1,"AGE",], mean = l.inputs[["AGE_start_mean"]], sd = l.inputs[["AGE_start_sd"]]),0)
  a.out[1,"AGE",][a.out[1,"AGE",]<50] <- 50
  a.out[1,"AGE",][a.out[1,"AGE",]>99] <- 99
  a.out[1,"SEX",]       <- f.qcat(p_rand = a.random[1,"SEX",], p_cat = c(l.inputs[["p.SEX_start_male"]], 
                                                                         l.inputs[["p.SEX_start_female"]]), values = l.inputs[["v.SEX_val"]])
  a.out[1,"EDU",]       <- f.qcat(p_rand = a.random[1,"EDU",], p_cat = l.inputs[["p.EDU_start"]], values = l.inputs[["v.EDU_val"]])
  a.out[1,"RACEETH",]   <- f.qcat(p_rand = a.random[1,"RACEETH",], p_cat = l.inputs[["p.RACEETH_start"]], values = l.inputs[["v.RACEETH_val"]])
  a.out[1,"INCOME",]    <- f.qcat(p_rand = a.random[1,"INCOME",], p_cat = l.inputs[["p.INCOME_start"]], values = l.inputs[["v.INCOME_val"]])
  a.out[1,"MEDBUR",]    <- qbeta(p = a.random[1,"MEDBUR",], shape1 = 1.5, shape2 = 6) * max(l.inputs[["v.MEDBUR_val"]])
  a.out[1,"APOE4",]     <- f.qcat(p_rand = a.random[1,"APOE4",], p_cat = l.inputs[["p.APOE4_start"]], values = l.inputs[[".APOE4_val"]])
  
  
  a.out[1,"TX",]         <- 0
  a.out[1,"SYN",]        <- 0    # Everyone starts healthy
  a.out[1,"BHA",]        <- 0    # Will be assigned as people are tested
  a.out[1,"CDR_track",]  <- NA   # Will be determined based on attributes
  a.out[1,"CDRfast_sd1",] <- qnorm(p = a.random[1,"CDRfast_sd1",], mean = 0, sd = l.inputs[["r.CDRfast_sd1"]]) 
  #%>%
  # rescale(to = c(-0.9*l.inputs[["r.CDRfast_mean"]], 0.9*l.inputs[["r.CDRfast_mean"]]))
  a.out[1,"CDRslow_sd1",] <- qnorm(p = a.random[1,"CDRslow_sd1",], mean = 0, sd = l.inputs[["r.CDRslow_sd1"]]) 
  #%>%
  # rescale(to = c(-0.9*l.inputs[["r.CDRslow_mean"]], 0.9*l.inputs[["r.CDRslow_mean"]]))
  a.out[1,"CDR",]       <- 0    # CDR-SB true score -- everyone starts cognitively normal
  a.out[1,"CDR_obs",]   <- -9   # CDR-SB observed score -- will be assigned as people are tested
  a.out[1,"SEV",]       <- NA   # Will be assigned as people are tested
  a.out[1,"SEV_obs",]   <- -9   # Will be assigned as people are tested
  a.out[1,"FUN",]       <- 0    # FAQ -- everyone starts with no functional impairment
  a.out[1,"BEH",]       <- NA
  a.out[1,"INSTIT",]    <- 0
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
      hr.mort_age   = l.inputs[["hr.mort_age"]]
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
      v.MEDBUR.lag    = a.out[t-1,"MEDBUR",alive],
      v.INCOME.lag    = a.out[t-1,"INCOME",alive],
      v.SYN.lag       = a.out[t-1,"SYN",alive], 
      random_cycle    = a.random[t,"SYN",alive], 
      n.alive         = n.alive
    )
    
    # CDR_track
    a.out[t,"CDR_track",alive] <- f.update_CDR_track(
      v.SEV.lag        = a.out[t-1,"SEV",alive],
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
      v.CDR.lag        = a.out[t-1,"CDR",alive], 
      r.CDRfast_mean   = l.inputs[["r.CDRfast_mean"]],
      r.CDRslow_mean   = l.inputs[["r.CDRslow_mean"]],
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
    
    # BHA
    a.out[t,"BHA",alive] <- f.update_BHA(
      v.SYN            = a.out[t,"SYN",alive],
      v.BHA.lag        = a.out[t-1,"BHA",alive],
      v.SEV            = a.out[t,"SEV",alive],
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
  a.out[a.out == -9] <- NA
  
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
  l.out[["healthy.trc"]] <- as.matrix(apply(X = a.out[,"SYN",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
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
  
  # state reside time by age at onset
  onset_cycle <- apply(a.out[,"SEV",], 2, function(x) {
    # Identify indices where SEV == 0 (MCI), ignoring NA
    mci_indices <- which(!is.na(x) & x == 0)
    
    # Return the first such index, or NA if none exist
    if (length(mci_indices) > 0) mci_indices[1] else NA
  })
  age_at_mci <- ifelse(!is.na(onset_cycle), a.out[cbind(onset_cycle, match("AGE", dimnames(a.out)[[2]]), seq_len(dim(a.out)[3]))], NA)
  age_bins <- cut(age_at_mci, breaks = seq(50, 100, by = 5), right = FALSE, include.lowest = TRUE)
  age_groups <- levels(age_bins)
  severity_levels <- c("mci", "mil", "mod", "sev")
  result_matrix <- matrix(NA, nrow = length(age_groups), ncol = length(severity_levels),
                          dimnames = list(age_groups, severity_levels))
  for (sev in 0:3) {
    # Total time in the current state for each individual
    time_in_state <- colSums(a.out[,"SEV",] == sev, na.rm = TRUE)
    
    # Average time per age group
    result_matrix[, sev + 1] <- tapply(time_in_state, age_bins, function(x) round(mean(x, na.rm = TRUE), digits = 2))
  }
  result_df <- as.data.frame(result_matrix) %>% 
    mutate(age_group = factor(age_groups,
                              labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+"))) %>%
    select(age_group, mci, mil, mod, sev) 
  rownames(result_df) <- NULL
  l.out[["reside_time_by_age_of_onset"]] <- result_df
  
  
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
  l.out[["state_concordance"]] <- matrix(data = NA, nrow = l.inputs[["n.cycle"]], ncol = 10,
                                         dimnames = list(NULL, c("h_h","h_mci","h_dem","mci_h","mci_mci","mci_dem","dem_h","dem_mci","dem_dem","dth")))
  l.out[["state_concordance"]][,"h_h"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"h_mci"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"h_dem"] <- as.matrix(apply(X = a.out[,"SYN",]==0 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]>=1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_h"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_mci"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"mci_dem"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]==0 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]>=1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_h"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_mci"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["state_concordance"]][,"dem_dem"] <- as.matrix(apply(X = a.out[,"SYN",]==1 & a.out[,"SEV",]>=1 & a.out[,"BHA",]==1 & a.out[,"SEV_obs",]>=1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
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
      at_risk <- which(a.out[t - 1, "AGE", ] >= age_min & a.out[t - 1, "AGE", ] <= age_max & a.out[t - 1, "SYN", ] == 0)
      
      # Calculate person-years at risk for this time point
      person_years <- person_years + length(at_risk)
      
      # Count new MCI cases (healthy in previous cycle and MCI in current cycle)
      new_cases <- new_cases + sum(a.out[t, "SYN", at_risk] == 1, na.rm = TRUE)
    }
    
    incidence_rate <- if (person_years > 0) new_cases / person_years else NA
    l.out[["mci_incidence"]][i, 1] <- incidence_rate
  }
  
  # Prevalence by severity by age -- this requires that age and cycle are equal (everyone starts the same age)
  l.out[["prevalence_grouped"]] <- as.data.frame(l.out[["state_trace"]]) %>%
    mutate(age = 50:99) %>%
    mutate(alive = 1 - dth,
           age_group = factor(cut(age, breaks = seq(50, 100, by = 5), right = FALSE),
                              labels = c("50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+"))) %>%
    mutate(mci = mci / alive,
           mil = mil / alive,
           mod = mod / alive,
           sev = sev / alive) %>%
    group_by(age_group) %>%
    summarize(avg_mci = mean(mci),
              avg_mil = mean(mil),
              avg_mod = mean(mod),
              avg_sev = mean(sev))
  
  
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
f.out_summary <- function(l.inputs, printLevel = 0) {
  
  # strategy 1
  l.inputs_strat1 <- l.inputs
  l.inputs_strat1[["strategy"]] <- l.inputs[["strategy_strat1"]]
  l.inputs_strat1[["Tx"]] <- l.inputs[["Tx_strat1"]]
  a.out_strat1 <- f.run(l.inputs = l.inputs_strat1, printLevel = printLevel)
  a.out_qc_strat1 <- f.qaly_cost(a.out = a.out_strat1, l.inputs = l.inputs_strat1)
  out_strat1 <- f.out_aggregate(a.out = a.out_qc_strat1, l.inputs = l.inputs_strat1)
  
  # strategy 2
  l.inputs_strat2 <- l.inputs
  l.inputs_strat2[["strategy"]] <- l.inputs[["strategy_strat2"]]
  l.inputs_strat2[["Tx"]] <- l.inputs[["Tx_strat2"]]
  a.out_strat2 <- f.run(l.inputs = l.inputs_strat2, printLevel = printLevel)
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
    l.inputs = l.inputs,
    a.out_strat1 = a.out_strat1, a.out_strat2 = a.out_strat2,
    a.out_qc_strat1 = a.out_qc_strat1, a.out_qc_strat2 = a.out_qc_strat2,
    out_strat1 = out_strat1, out_strat2 = out_strat2,
    m.out = m.out
  ))
  
}

# generate figures
f.figures <- function(l.out, l.inputs) {
  
  # Reshape syndrome data to long format for use with ggplot.
  
  #### True (unobserved) status ####
  t.trace_syndrome_true <- as.data.frame(l.out$state_trace[ ,c("healthy", "mci", "mil", "mod", "sev", "dth")],) %>%
    mutate(year = 1:l.inputs[["n.cycle"]] + 49)
  t.trace_syndrome_true_long <- t.trace_syndrome_true %>%
    pivot_longer(cols = -year, names_to = "syndrome", values_to = "proportion") %>%
    mutate(syndrome = fct_rev(factor(syndrome, levels = c("healthy", "mci", "mil", "mod", "sev", "dth"))))
  
  fig.progression_true <- ggplot(t.trace_syndrome_true_long, aes(x = year, y = proportion, fill = syndrome)) +
    geom_area(alpha = 0.8, position = "stack") + 
    geom_path(aes(group = syndrome), position = "stack", color = "black", linewidth = 0.5) + 
    scale_x_continuous(breaks =  seq(min(t.trace_syndrome_true_long$year), max(t.trace_syndrome_true_long$year), by = 5), minor_breaks = NULL) +
    scale_fill_manual(values = c("white","orange","purple","green","yellow","pink"),
                      labels = c("Death","Severe dementia","Moderate dementia","Mild dementia","MCI","Cognitively intact"),
                      guide = guide_legend(override.aes = list(colour = "black", size = 0.5))) +
    labs(title = "GRAM: Progression of Cognitive Impairment",
         subtitle = paste0(l.inputs[["scenario"]], "\n(N = ", l.inputs[["n.ind"]], " individuals)"),
         x = "Age",
         y = "Proportion of Population",
         fill = "True Cognitive Status") +
    theme(axis.text = element_text(size = 14),  # Increase axis tick font size
          axis.title = element_text(size = 16),
          panel.grid.major = element_line(color = "gray50", linewidth = 0.8),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  #### Observed status ####
  t.trace_syndrome_obs <- as.data.frame(l.out$state_concordance) %>%
    mutate(year = 1:l.inputs[["n.cycle"]] + 49)
  t.trace_syndrome_obs_long <- t.trace_syndrome_obs %>%
    pivot_longer(cols = -year, names_to = "status", values_to = "proportion") %>%
    mutate(status = fct_rev(factor(status, levels = c("h_h","mci_h","dem_h","h_mci","mci_mci","dem_mci","h_dem","mci_dem","dem_dem","dth"))))
  
  t.trace_syndrome_obs_long <- t.trace_syndrome_obs_long %>%
    mutate(
      pattern = case_when(
        status %in% c("h_h","mci_mci","dem_dem","dth") ~ "none",
        status %in% c("mci_h","dem_h","dem_mci") ~ "stripe",
        status %in% c("h_mci","mci_dem","h_dem") ~ "circle"
      ),
      observed_color = case_when(
        status %in% c("h_h","mci_h","dem_h") ~ "pink",
        status %in% c("h_mci","mci_mci","dem_mci") ~ "yellow",
        status %in% c("h_dem","mci_dem","dem_dem") ~ "darkgreen",
        status == "dth" ~ "white"
      )
    )
  
  
  fig.progression_obs <- ggplot(t.trace_syndrome_obs_long, aes(x = year, y = proportion, fill = observed_color, pattern = pattern)) +
    #  geom_area(alpha = 0.8, position = "stack") + 
    #  geom_path(aes(group = status), position = "stack", color = "black", linewidth = 0.5) + 
    geom_area_pattern(color = "black",
                      alpha = 0.8,
                      pattern_fill = "black",
                      pattern_density = 0.1,
                      pattern_spacing = 0.01) +
    scale_x_continuous(breaks =  seq(min(t.trace_syndrome_obs_long$year), max(t.trace_syndrome_obs_long$year), by = 5), minor_breaks = NULL) +
    scale_fill_identity(name = "Observed Status",
                        labels = c("dth","dem","mci","h"),
                        guide = guide_legend(override.aes = list(colour = "black", size = 0.5))) +
    scale_pattern_manual(values = c("none" = "none", "circle" = "circle", "stripe" = "stripe"),
                         name = "Concordance") +
    labs(title = "GRAM: Progression of Cognitive Impairment",
         subtitle = paste0(l.inputs[["scenario"]], "\n(N = ", l.inputs[["n.ind"]], " individuals)"),
         x = "Age",
         y = "Proportion of Population",
         fill = "True Status",
         pattern = "Observed Status") +
    theme(axis.text = element_text(size = 14),  # Increase axis tick font size
          axis.title = element_text(size = 16),
          panel.grid.major = element_line(color = "gray50", linewidth = 0.8),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          title = element_text(size = 16))
  
  
  return(list(
    fig.progression_true = fig.progression_true,
    fig.progression_obs = fig.progression_obs))
}

######################################## 3.3. MAIN FUNCTIONS: WRAPPER  ########################################

f.wrap_run <- function(l.inputs, printLevel = 0) {
  
  
  output <- f.run(l.inputs = l.inputs, printLevel = printLevel)
  
  # Aggregate results for full cohort
  aggregated_results_totpop <- f.out_aggregate(a.out = output, l.inputs = l.inputs)
  fig.progression <- f.figures(l.out = aggregated_results_totpop, l.inputs = l.inputs)
  
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

