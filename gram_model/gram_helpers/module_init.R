######################################## GRAM MODULE: INITIALIZATION ########################################
# These functions initialize the model. 

f.initialize <- function(l.inputs, microdata) {
  
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
    
    synthetic_pop <- generate_synthetic_sample(microdata, target_size = l.inputs[["n.ind"]], weights = microdata$PERWT, seed = l.inputs[["seed_stochastic"]])
    
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
  a.out[1,"CDR",]         <- case_match(a.out[1,"SYN",],
                                        0 ~ 0,
                                        1 ~ qunif(p = a.random[1,"CDR",], min = l.inputs[["cutoff_CDR"]]["mci"], max = l.inputs[["cutoff_CDR"]]["moderate"]))
  a.out[1,"MEMLOSS",]   <- case_match(a.out[1,"SYN",],
                                      0 ~ 0,
                                      1 ~ f.qcat(p_rand = a.random[1,"MEMLOSS",], p_cat = l.inputs[["p.MEMLOSS_start"]], values = l.inputs[["v.MEMLOSS_val"]]))   # !! when starting with prevalent pop, this will cause a problem as people with dem might get MEMLOSS of 1
  a.out[1,"SEV",]       <- f.update_SEV(v.SYN = a.out[1,"SYN",], v.CDR = a.out[1,"CDR",], 
                                        cutoff_CDR = l.inputs[["cutoff_CDR"]], n.alive = l.inputs[["n.ind"]])   # Everyone starts healthy
  a.out[1,"COGCON",]    <- 0    # assume no concerns at age 50 when healthy !! TODO: make this dynamic!
  a.out[1,"BHA",]       <- -9   # will be assigned as people are tested
  a.out[1,"last_BHA_age",]   <- NA   # tracking variable for implementing BHA scenarios
  a.out[1,"any_BHA_pos",]    <- FALSE
  a.out[1,"CDR_obs",]   <- -9   # DR-SB observed score -- will be assigned as people are tested
  a.out[1,"SEV_obs",]   <- -9   # will be assigned as people are tested
  
  a.out[1,"DX",]        <- f.qcat(p_rand = a.random[1,"DX",], p_cat = l.inputs[["p.DX_start"]], values = l.inputs[["v.DX_val"]])
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