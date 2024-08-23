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
  a.out[1,"SYN",]       <- 1
  a.out[1,"SEV",]       <- NA    # Will be assigned as people develop dementia
  a.out[1,"COG",]       <- 4     # CDR-SB -- everyone starts cognitively normal
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
  
  # time in mild cognition
  l.out[["SEV1.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV1.sum"]] <- sum(l.out[["SEV1.trc"]])
  l.out[["SEV1.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV1.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV1.sum.dis"]] <- sum(l.out[["SEV1.trc.dis"]])
  # time in moderate cognition
  l.out[["SEV2.trc"]] <- as.matrix(apply(X = a.out[,"SEV",]==2, MARGIN = 1, FUN = sum, na.rm = TRUE)/l.inputs[["n.ind"]])
  l.out[["SEV2.sum"]] <- sum(l.out[["SEV2.trc"]])
  l.out[["SEV2.trc.dis"]] <- as.matrix(f.discount(x = l.out[["SEV2.trc"]], discount_rate = l.inputs[["r.discount_QALY"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["SEV2.sum.dis"]] <- sum(l.out[["SEV2.trc.dis"]])
  # time in severe cognition
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
    mutate(year = (1:l.inputs[["n.cycle"]]) - 1)
  t.trace_syndrome_long <- t.trace_syndrome %>%
    pivot_longer(cols = -year, names_to = "syndrome", values_to = "proportion") %>%
    mutate(syndrome = factor(syndrome, levels = c("healthy", "mci", "mil", "mod", "sev", "dth")))
  
  fig.syndrome_trace <- ggplot(t.trace_syndrome_long, aes(x = year, y = proportion, fill = fct_rev(syndrome))) +
    geom_area(alpha = 0.8, position = "stack") + 
    scale_x_continuous(breaks = unique(t.trace_syndrome_long$year), minor_breaks = NULL) +
    theme_minimal()
  
  return(fig.syndrome_trace)
}

######################################## 2.4. FUNCTIONS: WRAPPERS ########################################

f.wrap_run <- function(a.random, l.inputs) {
  output <- f.run(a.random = a.random, l.inputs = l.inputs)
  tables <- f.out_aggregate(a.out = output, l.inputs = l.inputs)
  figures <- f.figures(l.out = tables, l.inputs = l.inputs)
  
  print(figures)
  invisible(return(list(
    output = output,
    tables = tables,
    figures = figures)
  ))
}
