######################################## GRAM MODULE: MEDICAL RECORD ########################################
# This section defines the functions for updating known/observed health status variables (cognitive tests, neuropsychiatric assessments etc)

#### Module Wrapper ####
f.module_medical_record <- function(l.inputs, a.out, t, a.random, alive, n.alive) {
  
  # COGCON
  a.out[t,"COGCON",alive] <- f.update_COGCON(
    scenario         = l.inputs[["scenario"]],
    v.AGE            = a.out[t,"AGE",alive],
    v.SYN            = a.out[t,"SYN",alive],
    v.SEV            = a.out[t,"SEV",alive],
    v.COGCON.lag     = a.out[t-1,"COGCON",alive],
    rr.cogcon_prior  = l.inputs[["rr.cogcon_prior"]],
    v.DX.lag         = a.out[t-1,"DX",alive],
    random_cycle     = a.random[t,"COGCON",alive],
    n.alive          = n.alive
  )
  
  
  # BHA
  a.out[t,"BHA",alive] <- f.update_BHA(
    scenario         = l.inputs[["scenario"]],
    cycle            = t,
    v.HCARE          = a.out[t,"HCARE",alive],
    v.DX.lag         = a.out[t-1,"DX",alive],
    v.COGCON         = a.out[t,"COGCON",alive],
    v.AGE            = a.out[t,"AGE",alive],
    v.last_BHA_age   = a.out[t-1,"last_BHA_age",alive],
    v.BHA.lag        = a.out[t-1,"BHA",alive],
    v.NP.lag         = a.out[t-1,"NP",alive], 
    v.PET.lag        = a.out[t-1,"PET",alive],
    v.SYN            = a.out[t,"SYN",alive],
    v.SEV            = a.out[t,"SEV",alive],
    v.MEMLOSS        = a.out[t,"MEMLOSS",alive], 
    random_cycle     = a.random[t,"BHA",alive],
    n.alive          = n.alive,
    v.any_BHA_pos    = a.out[t-1,"any_BHA_pos",alive],
    v.last_FP_age    = v.last_FP_age[alive],
    v.any_PCP_pos    = a.out[t-1,"any_PCP_pos",alive]
  )
  
  a.out[t,"last_BHA_age",alive] <- ifelse(a.out[t,"BHA",alive] >= 0, a.out[t,"AGE", alive], a.out[t-1,"last_BHA_age", alive])
  a.out[t,"any_BHA_pos",alive] <-  as.numeric((a.out[t-1,"any_BHA_pos",alive]) | (a.out[t,"BHA",alive] == 1))
  
  
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
    v.CDR_obs      = a.out[t,"CDR_obs",alive],
    cutoff_CDR     = l.inputs[["cutoff_CDR"]],
    n.alive        = n.alive
  ) 
  
  
  # PCP
  a.out[t,"PCP",alive] <- f.update_PCP(
    prob_pcpfu         = l.inputs[["scenario"]][["prob_pcpfu"]],
    v.DX.lag           = a.out[t-1,"DX",alive],
    v.HCARE            = a.out[t,"HCARE",alive],
    v.BHA              = a.out[t,"BHA",alive],
    v.SYN              = a.out[t,"SYN",alive],
    v.SEV              = a.out[t,"SEV",alive],
    p.PCP_confirm_TP   = l.inputs[["scenario"]][["p.PCP_confirm_TP"]] %||% l.inputs[["p.PCP_confirm_TP"]],
    p.PCP_reject_FP    = l.inputs[["scenario"]][["p.PCP_reject_FP"]] %||% l.inputs[["p.PCP_reject_FP"]],
    random_cycle       = a.random[t,"PCP",alive],
    n.alive            = n.alive
  )
  
  a.out[t,"any_PCP_pos",alive] <-  as.numeric((a.out[t-1,"any_PCP_pos",alive]) | (a.out[t,"PCP",alive] == 1))
  
  
  # PET 
  a.out[t,"PET",alive] <- f.update_PET(
    v.PET.lag     = a.out[t-1,"PET",alive]
  )
  
  
  # NP
  a.out[t,"NP",alive] <- f.update_NP(
    v.NP.lag     = a.out[t-1,"NP",alive]
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
  
  return(a.out)
}



#### Module Functions ####
######################################## COGCON
f.update_COGCON <- function(scenario, v.AGE, v.SYN, v.SEV, v.COGCON.lag, rr.cogcon_prior, v.DX.lag, random_cycle, n.alive) {
  
  cogcon <- rep(-9, n.alive)
  
  if(!is.null(scenario[["test"]])) {
    
    select_col <- case_when(
      v.SYN < 1 ~ 2,
      v.SEV == 0 ~ 3,
      v.SEV >= 1 ~ 4
    )
    
    cogcon_lookup_coordinates <- matrix(data = c(round(v.AGE,0)-50+1, select_col), ncol = 2)  
    
    prob_cogcon <- scenario[["probs_cogcon"]][cogcon_lookup_coordinates]
    prob_cogcon[v.COGCON.lag == 1] <- f.adjustprobability(prob_cogcon[v.COGCON.lag == 1], t_new = 1, t_old = 1, RR = rr.cogcon_prior)
    
    cogcon[v.DX.lag == 0] <- as.numeric(prob_cogcon[v.DX.lag == 0] > random_cycle[v.DX.lag == 0])
  } 
  
  return(cogcon)
}



######################################## BHA


f.update_BHA <- function(scenario, cycle, v.HCARE, v.DX.lag, v.COGCON, v.AGE, v.last_BHA_age, v.BHA.lag, 
                         v.NP.lag = NULL, v.PET.lag = NULL, v.SYN, v.SEV, v.MEMLOSS, random_cycle, n.alive, 
                         v.any_BHA_pos, v.last_FP_age, v.any_PCP_pos) {
  
  bha <- rep(-9, n.alive)
  
  if(!is.null(scenario[["test"]])) {
    
    # universally eligible for assessment
    assess <- (v.HCARE == scenario$HCARE) & (v.DX.lag == scenario$DX)
    
    # age criteria for assessment
    assess <- assess & (v.AGE >= scenario$age_first_test)
    
    # cohort split criteria for initial assessment
    if(is.null(scenario$cohort_split)) {scenario$cohort_split <- 1}
    first_test <- is.na(v.last_BHA_age) & (((as.numeric(names(v.last_BHA_age))-1) %% scenario$cohort_split) + 1) == (((cycle-1) %% scenario$cohort_split) + 1)
    
    # repeat criteria for assessment
    interval_ok <- !is.na(v.last_BHA_age) & ((v.AGE - v.last_BHA_age) >= scenario$repeat_interval)
    assess <- assess & (first_test | interval_ok)
    
    # optional stop rule
    stop_test <- scenario$stop_rule(any_BHA_pos = v.any_BHA_pos, 
                                    NP = v.NP.lag, 
                                    PET = v.PET.lag,
                                    any_PCP_pos = v.any_PCP_pos, 
                                    repeat_after_FP = scenario$repeat_after_FP)
    assess <- assess & (is.na(stop_test) | !stop_test) & (v.AGE <= scenario$age_stop_test)
    
    # cognitive concerns criteria
    eligible <- assess & (v.COGCON == 1)
    
    # neuropsych criteria
    if (!is.null(scenario$NP)) {
      eligible <- eligible & (is.null(v.NP.lag) | v.NP.lag != 1)
    }
    
    # pause after detected false positive BHA
    if (!is.null(scenario$repeat_after_FP)) {
      interval_after_fp_ok <- is.na(v.last_FP_age) | ((v.AGE - v.last_FP_age) >= scenario$repeat_after_FP)
      eligible <- eligible & interval_after_fp_ok
    }
    
    # calculate test results
    bha[assess & !eligible] <- -8
    
    if(any(eligible)) {
      bha[eligible] <- case_when(
        v.SYN[eligible] == 0      ~ as.numeric((1 - scenario$specificity) > random_cycle[eligible]),
        v.SYN[eligible] == 0.5 & (v.BHA.lag[eligible] >= 0)  ~ as.numeric((scenario$sensitivity[1] * 1) > random_cycle[eligible]),  # TCI and second consecutive BHA (sens is lower due to practice effect)
        v.SYN[eligible] == 0.5 & (v.BHA.lag[eligible] < 0)  ~ as.numeric(scenario$sensitivity[1] > random_cycle[eligible]),  # TCI and first BHA
        v.MEMLOSS[eligible] == 1  ~ as.numeric(scenario$sensitivity[2] > random_cycle[eligible]),
        v.SEV[eligible] == 0      ~ as.numeric(scenario$sensitivity[3] > random_cycle[eligible]),
        v.SEV[eligible] >= 1      ~ as.numeric(scenario$sensitivity[4] > random_cycle[eligible])
      )
    }
  }
  
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

######################################## DX - outside of the BHA pathway

f.update_DX <- function(v.DX.lag, v.SYN, v.SEV, v.HCARE, random_cycle, n.alive) {
  dx <- rep(0, n.alive)
  
  # if already diagnosed, stays diagnosed
  dx[v.DX.lag == 1] <- 1 
  
  # if previously undiagnosed, check SEV in last cycle
  mci <- v.SYN == 1 & v.SEV == 0 & !is.na(v.SEV) & v.HCARE == 1
  mil <- v.SYN == 1 & v.SEV == 1 & !is.na(v.SEV) & v.HCARE == 1 
  mod_sev <- v.SYN == 1 & v.SEV %in% c(2,3) & !is.na(v.SEV) & v.HCARE == 1
  
  dx[v.DX.lag == 0 & mci] <- as.numeric(random_cycle[v.DX.lag == 0 & mci] < 0.15)
  dx[v.DX.lag == 0 & mil] <- as.numeric(random_cycle[v.DX.lag == 0 & mil] < 0.40)
  dx[v.DX.lag == 0 & mod_sev] <- as.numeric(random_cycle[v.DX.lag == 0 & mod_sev] < 0.85)
  
  return(dx)
}


######################################## PCP

f.update_PCP <- function(prob_pcpfu, v.DX.lag, v.HCARE, v.BHA, v.SYN, v.SEV,
                         p.PCP_confirm_TP, p.PCP_reject_FP, random_cycle, n.alive) {

  pcp <- rep(NA, n.alive)

  if(is.null(prob_pcpfu) || prob_pcpfu == 0) {
    eligible <- rep(0, n.alive)
  } else {
    eligible <- v.DX.lag == 0 & v.HCARE == 1 & v.BHA == 1
    select_fu <- rbinom(n = length(eligible), size = 1, prob = prob_pcpfu)
    eligible <- eligible & select_fu
  }

  # PCP as filter on BHA-positive patients.
  # Truly healthy: PCP correctly rejects FP BHA with probability p_PCP_reject_FP
  pcp[eligible & v.SYN < 1] <- as.numeric((1 - p.PCP_reject_FP) > random_cycle[eligible & v.SYN < 1])

  # Truly impaired: PCP confirms TP BHA with probability varying by SEV (SEV 0-3 -> index 1-4)
  if(any(eligible & v.SYN == 1)) {
    confirm_prob <- p.PCP_confirm_TP[v.SEV[eligible & v.SYN == 1] + 1]
    pcp[eligible & v.SYN == 1] <- as.numeric(confirm_prob > random_cycle[eligible & v.SYN == 1])
  }

  return(pcp)
}


######################################## PET

f.update_PET <- function(v.PET.lag) {
  pet <- v.PET.lag
  return(pet)
}


######################################## NP

f.update_NP <- function(v.NP.lag) {
  np <- v.NP.lag
  return(np)
}


