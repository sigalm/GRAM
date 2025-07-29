######################################## GRAM MAIN FUNCTIONS: RUN MODEL ########################################

# run model
f.run <- function(l.inputs, microdata, printLevel) {
  
  # progress bar
  if (printLevel > 1) {
    ptm <- proc.time()
    stime <- Sys.time()
    pb = txtProgressBar(min = 0, max = l.inputs[["n.cycle"]], initial = 0)
  }
  
  init <- f.initialize(l.inputs = l.inputs, microdata = microdata)
  a.random <- init$a.random
  a.out <- init$a.out
  
  
  # update progress bar
  if (printLevel > 1) setTxtProgressBar(pb, 1)
  
  # run subsequent cycles
  for(t in 2:l.inputs[["n.cycle"]]) {
    
    # TIME
    a.out[t,"TIME",] <- a.out[t-1,"TIME",] + 1

    ########## !!!!!!!!!!!!!!! The model is run by updating each attribute in a loop over the cycles (over time). 
    # At each cycle attributes are updated using the attribute status at the previous cycle or the status at the current cycle. 
    # Except the first cycle, which is manually put in (i.e., starting values). 
    # For transparency, no information from other than the previous or current is used. Information from more than 1 cycle ago 
    #     could be used by tracking the history of an attribute in a separate attribute. For each attribute a function is written to update it. 
    #     Then, in a loop all functions are called to update their status cycle by cycle. 
    
    a.out <- f.module_mortality(a.out, t, a.random)
    
    # identify those alive at current observation (to be used for subsetting the other functions 
    #       so they don't have to process the data of the individuals no longer alive)
    alive <- a.out[t,"ALIVE",]==1
    n.alive <- sum(alive)   # number of individuals alive
    
    a.out <- f.module_socdem(a.out, t, a.random, alive)

    a.out <- f.module_true_health(a.out, t, a.random, alive, n.alive)
    
    a.out <- f.module_medical_record(a.out, t, a.random, alive, n.alive)

    a.out <- f.module_treatment(a.out, t, a.random, alive, n.alive)
      
    # update progress bar
    if (printLevel > 1) setTxtProgressBar(pb, t)
    
  }
  
  # The "round(score*2)/2" ensures that scores are increments of 0.5 
  a.out[,"CDR",] <- ifelse(a.out[,"CDR",] >= 0.5, (round(a.out[,"CDR",] * 2) / 2), 0)
  a.out[,"CDR_obs",] <- ifelse(a.out[,"CDR",] >= 0.5, (round(a.out[,"CDR_obs",] * 2) / 2), 0)
  
  # run time
  if (printLevel >1 ) {
    print(proc.time() - ptm)
    print(Sys.time() - stime)
  }
  
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
  
  # COST: BHA
  COST_test <- as.numeric(a.out[,"BHA",]!=-9) * l.inputs[["c.bha"]] + as.numeric(a.out[,"BHA",] == 1) * l.inputs[["c.bhapos"]]
  COST_test[is.na(COST_test)] <- 0
  
  # COST: PET
  COST_pet <- as.numeric(a.out[,"PET",]) * l.inputs[["c.pet"]]
  COST_pet[is.na(COST_pet)] <- 0
  
  # COST: NP Assessment
  COST_np <- as.numeric(a.out[,"NP",]) * l.inputs[["c.np"]]
  COST_np[is.na(COST_np)] <- 0
  
  # COST: treatment (DMT)
  COST_tx <- as.numeric(a.out[,"TX",]) * l.inputs[["c.Tx"]]
  COST_tx[is.na(COST_tx)] <- 0
  
  # COST: treatment (non-DMT)
  COST_tx2 <- as.numeric(a.out[,"TX2",]) * l.inputs[["c.Tx2"]]
  COST_tx2[is.na(COST_tx2)] <- 0
  
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
  a.out[,"COST_test",] <- COST_test
  a.out[,"COST_fu",] <- COST_pet + COST_np
  a.out[,"COST_tx",] <- COST_tx
  a.out[,"COST_tx2"] <- COST_tx2
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
  l.out[["COST_test"]] <- as.matrix(apply(X = a.out[,"COST_test",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["COST_test.sum"]] <- sum(l.out[["COST_test"]])
  l.out[["COST_test.dis"]] <- as.matrix(f.discount(x = l.out[["COST_test"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_test.dis.sum"]] <- sum(l.out[["COST_test.dis"]])
  
  # COST_tx
  l.out[["COST_tx"]] <- as.matrix(apply(X = a.out[,"COST_tx",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["COST_tx.sum"]] <- sum(l.out[["COST_tx"]])
  l.out[["COST_tx.dis"]] <- as.matrix(f.discount(x = l.out[["COST_tx"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_tx.dis.sum"]] <- sum(l.out[["COST_tx.dis"]])
  
  # COST_fu
  l.out[["COST_fu"]] <- as.matrix(apply(X = a.out[,"COST_fu",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n)
  l.out[["COST_fu.sum"]] <- sum(l.out[["COST_fu"]])
  l.out[["COST_fu.dis"]] <- as.matrix(f.discount(x = l.out[["COST_fu"]], discount_rate = l.inputs[["r.discount_COST"]], n.cycle = l.inputs[["n.cycle"]]))
  l.out[["COST_fu.dis.sum"]] <- sum(l.out[["COST_fu.dis"]])
  
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

