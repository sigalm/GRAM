######################################## GRAM MODULE: MEDICAL RECORD ########################################
# This section defines the functions for updating known/observed health status variables (cognitive tests, neuropsychiatric assessments etc)

#### Module Wrapper ####
f.module_medical_record <- function(l.inputs, a.out, t, a.random, alive, n.alive) {
  
  # Who CAN be tested this cycle: has a provider, is undiagnosed, is in the age window,
  # is due under the cohort split / repeat interval, and has not hit a stop condition.
  # Computed once here because selection is drawn among these people, and only these.
  assess <- f.assess_eligible(
    scenario         = l.inputs[["scenario"]],
    cycle            = t,
    v.HCARE          = a.out[t,"HCARE",alive],
    v.DX.lag         = a.out[t-1,"DX",alive],
    v.AGE            = a.out[t,"AGE",alive],
    v.last_BHA_age   = a.out[t-1,"last_BHA_age",alive],
    v.NP.lag         = a.out[t-1,"NP",alive],
    v.any_BHA_pos    = a.out[t-1,"any_BHA_pos",alive],
    v.any_PCP_pos    = a.out[t-1,"any_PCP_pos",alive],
    n.alive          = n.alive
  )

  # SELECT
  a.out[t,"SELECT",alive] <- f.update_SELECT(
    scenario         = l.inputs[["scenario"]],
    assess           = assess,
    v.AGE            = a.out[t,"AGE",alive],
    v.SYN            = a.out[t,"SYN",alive],
    v.SEV            = a.out[t,"SEV",alive],
    v.SELECT.lag     = a.out[t-1,"SELECT",alive],
    rr.select_prior  = l.inputs[["scenario"]][["rr.select_prior"]],
    v.BHA.lag        = a.out[t-1,"BHA",alive],
    v.DX.lag         = a.out[t-1,"DX",alive],
    random_cycle     = a.random[t,"SELECT",alive],
    n.alive          = n.alive
  )
  
  
  # BHA
  a.out[t,"BHA",alive] <- f.update_BHA(
    scenario         = l.inputs[["scenario"]],
    assess           = assess,
    v.SELECT         = a.out[t,"SELECT",alive],
    v.declined_last  = f.declined_last_offer(a.out, t, alive),
    v.AGE            = a.out[t,"AGE",alive],
    v.BHA.lag        = a.out[t-1,"BHA",alive],
    v.SYN            = a.out[t,"SYN",alive],
    v.SEV            = a.out[t,"SEV",alive],
    v.MEMLOSS        = a.out[t,"MEMLOSS",alive], 
    random_accept    = a.random[t,"ACCEPT",alive],
    random_cycle     = a.random[t,"BHA",alive],
    n.alive          = n.alive
  )
  
  # last_BHA_age is the age at the last OFFER: a test taken (BHA 0/1) or one declined
  # (selected, due, and still not tested, BHA -8 with SELECT 1). repeat_interval counts
  # from it, so someone who declines is not offered again until the interval has passed,
  # exactly as if they had tested. Where probs_accept is unset nobody declines, and this
  # is the age at the last test, as it always was.
  offered <- (a.out[t,"BHA",alive] >= 0) | (a.out[t,"BHA",alive] == -8 & a.out[t,"SELECT",alive] == 1)
  a.out[t,"last_BHA_age",alive] <- ifelse(offered, a.out[t,"AGE", alive], a.out[t-1,"last_BHA_age", alive])
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
    p.PCP_confirm_TP   = l.inputs[["scenario"]][["p.PCP_confirm_TP"]],
    p.PCP_reject_FP    = l.inputs[["scenario"]][["p.PCP_reject_FP"]],
    random_cycle       = a.random[t,"PCP",alive],
    n.alive            = n.alive
  )
  
  a.out[t,"any_PCP_pos",alive] <-  as.numeric((a.out[t-1,"any_PCP_pos",alive]) | (a.out[t,"PCP",alive] == 1))
  
  
  
  
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
######################################## SELECT
f.update_SELECT <- function(scenario, assess, v.AGE, v.SYN, v.SEV, v.SELECT.lag,
                            rr.select_prior, v.BHA.lag, v.DX.lag, random_cycle, n.alive) {

  selected <- rep(-9, n.alive)

  if(!is.null(scenario[["test"]])) {

    # SELECT records selection status AS OF THE LAST ASSESSMENT, so it is carried forward
    # through cycles in which no assessment could take place. Redrawing every cycle would
    # accumulate selections in years the programme was never going to test anyone. The
    # carry-forward is also what select_persists reads to know who is already flagged.
    selected <- v.SELECT.lag

    draw <- assess & (v.DX.lag == 0)

    # An EHR-derived flag is a property of the record, not a judgement the patient
    # remakes every year: once the algorithm has flagged someone it goes on flagging
    # them, whatever the last test showed. Only the already-flagged are held out of the
    # redraw -- anyone not flagged draws at probs_select at every assessment, exactly as
    # they otherwise would, so this raises the FLOOR without touching the entry rate.
    # Behavioural selection -- a concern raised, a question endorsed, an opt-in -- is
    # redrawn in full at every assessment, which is the default.
    if (isTRUE(scenario[["select_persists"]])) draw <- draw & (v.SELECT.lag != 1)

    if (any(draw)) {
      prob_select <- f.lookup_by_state(scenario[["probs_select"]], v.AGE, v.SYN, v.SEV)

      # Optional: adjust the probability by the RESULT of a test in the immediately
      # preceding cycle, not merely by having been selected. A negative result reassures,
      # so the person is less likely to raise a concern (or endorse one) again; a positive
      # result is a different situation and carries its own value.
      #
      # Keyed on last cycle's BHA rather than on a carried-forward last result because
      # "no test in the prior cycle" is deliberately a redraw at the unadjusted
      # probability: the effect is immediate, and is not meant to persist across the
      # years between assessments. In the 3-yearly arms no one is ever assessed in
      # consecutive cycles, so the adjustment is inert there, which is the intended
      # reading -- reassurance does not survive a 3-year gap. BHA.lag is -9 (not
      # assessed) or -8 (assessed but not tested) when there is no result; both leave
      # the probability alone.
      rr <- f.rr_by_prior_result(rr.select_prior, v.BHA.lag)
      bump <- draw & (rr != 1)
      if (any(bump)) {
        prob_select[bump] <- f.adjustprobability(prob_select[bump], t_new = 1, t_old = 1, RR = rr[bump])
      }

      selected[draw] <- as.numeric(prob_select[draw] > random_cycle[draw])
    }

    # A prior diagnosis ends selection permanently (DX never reverts).
    selected[v.DX.lag == 1] <- -9
  }

  return(selected)
}



######################################## ASSESSMENT ELIGIBILITY

# Who CAN be tested this cycle, before selection is considered. Lifted out of
# f.update_BHA so that f.update_SELECT can draw among exactly this group: selection is
# a property of an assessment, so it should only be drawn when an assessment happens.
f.assess_eligible <- function(scenario, cycle, v.HCARE, v.DX.lag, v.AGE, v.last_BHA_age,
                              v.NP.lag, v.any_BHA_pos, v.any_PCP_pos, n.alive) {

  if (is.null(scenario[["test"]])) return(rep(FALSE, n.alive))

  # The no-prior-diagnosis condition is fixed, not a scenario choice.
  assess <- (v.HCARE == scenario$HCARE) & (v.DX.lag == 0)

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
                                  any_PCP_pos = v.any_PCP_pos)
  assess <- assess & (is.na(stop_test) | !stop_test) & (v.AGE <= scenario$age_stop_test)

  # neuropsych criteria. Checked here rather than in f.update_BHA so that everyone who is
  # selected and due is genuinely offered a test: SELECT 1 with BHA -8 then means a
  # decline and nothing else.
  if (!is.null(scenario$NP)) {
    assess <- assess & (is.null(v.NP.lag) | v.NP.lag != 1)
  }

  assess
}


# Did each living person decline the last test they were offered? Worked out from the
# history rather than carried in an attribute. last_BHA_age is the age at the last offer,
# and AGE rises by exactly one a cycle, so the offer was (AGE - last_BHA_age) cycles ago;
# a -8 there is a decline, since an offer is either taken (0/1) or declined (-8). Read
# from t-1's last_BHA_age, so this cycle's offer is never its own history.
f.declined_last_offer <- function(a.out, t, alive) {
  idx      <- which(alive)
  last_age <- a.out[t-1, "last_BHA_age", idx]
  declined <- rep(FALSE, length(idx))
  has      <- !is.na(last_age)
  if (any(has)) {
    cyc <- t - (a.out[t, "AGE", idx[has]] - last_age[has])
    bha_col <- match("BHA", dimnames(a.out)[[2]])
    declined[has] <- a.out[cbind(cyc, bha_col, idx[has])] == -8
  }
  declined
}


######################################## BHA


f.update_BHA <- function(scenario, assess, v.SELECT, v.declined_last, v.AGE, v.BHA.lag,
                         v.SYN, v.SEV, v.MEMLOSS, random_accept, random_cycle, n.alive) {
  
  bha <- rep(-9, n.alive)
  
  if(!is.null(scenario[["test"]])) {
    
    # Offered a test: due this cycle (assess) and selected (drawn in f.update_SELECT
    # among exactly this assess group).
    offered <- assess & (v.SELECT == 1)

    # Whether an offer is taken up. probs_accept is P(test | selected, true state), so
    # probs_select x probs_accept is P(tested) at a first offer. Someone who declined
    # the last offer they had takes the next one up at p.accept_after_decline instead,
    # whatever their state. Without probs_accept every offer is taken up and selection
    # alone decides who is tested, as it did before selection and acceptance were split.
    #
    # The acceptance draw has its own random slot. Sharing random_cycle with the test
    # result would tie the two together: only people with a low draw would accept, and
    # the same low draw would then make a positive result more likely.
    tested <- offered
    if (!is.null(scenario[["probs_accept"]]) && any(offered)) {
      prob_accept <- f.lookup_by_state(scenario[["probs_accept"]], v.AGE, v.SYN, v.SEV)
      if (!is.null(scenario[["p.accept_after_decline"]])) {
        prob_accept[v.declined_last] <- scenario[["p.accept_after_decline"]]
      }
      tested <- offered & (prob_accept > random_accept)
    }
    
    # Due but not tested: -8. SELECT says why -- 0 not selected, 1 selected and declined.
    bha[assess & !tested] <- -8
    
    if(any(tested)) {
      bha[tested] <- case_when(
        v.SYN[tested] == 0      ~ as.numeric((1 - scenario$specificity) > random_cycle[tested]),
        v.SYN[tested] == 0.5 & (v.BHA.lag[tested] >= 0)  ~ as.numeric((scenario$sensitivity[1] * 1) > random_cycle[tested]),  # TCI and second consecutive BHA (sens is lower due to practice effect)
        v.SYN[tested] == 0.5 & (v.BHA.lag[tested] < 0)  ~ as.numeric(scenario$sensitivity[1] > random_cycle[tested]),  # TCI and first BHA
        v.MEMLOSS[tested] == 1  ~ as.numeric(scenario$sensitivity[2] > random_cycle[tested]),
        v.SEV[tested] == 0      ~ as.numeric(scenario$sensitivity[3] > random_cycle[tested]),
        v.SEV[tested] >= 1      ~ as.numeric(scenario$sensitivity[4] > random_cycle[tested])
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
    if (is.null(p.PCP_confirm_TP) || is.null(p.PCP_reject_FP)) {
      stop("Scenario requests PCP follow-up (prob_pcpfu = ", prob_pcpfu,
           ") but p.PCP_confirm_TP / p.PCP_reject_FP are not set in the scenario config.",
           call. = FALSE)
    }
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


######################################## NP

f.update_NP <- function(v.NP.lag) {
  np <- v.NP.lag
  return(np)
}


