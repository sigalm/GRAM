######################################## GRAM MODULE: TREATMENT AND CARE ########################################
# This section defines the functions for updating treatment and dementia care variables.

#### Module Wrapper ####
f.module_treatment <- function(a.out, t, a.random, alive, n.alive) {
  
  # TX2
  a.out[t,"TX2",alive] <- f.update_TX2(
    v.TX2.lag     = a.out[t-1,"TX2",alive]
  )
  
  # LTC       # care costs calculated using average cost!
  a.out[t,"LTC",alive] <- f.update_LTC(
    v.LTC.lag     = a.out[t-1,"LTC",alive]
  )
  
  
  # TX
  

  a.out[t,"TX",alive] <- f.update_TX(
    v.SEV_obs            = a.out[t,"SEV_obs",alive], 
    v.SEV_obs.lag        = a.out[t-1,"SEV_obs",alive],
    v.TX.lag             = a.out[t-1,"TX",alive], 
    n.alive              = n.alive, 
    random_cycle         = a.random[t,"TX",alive], 
    Tx_t_max             = l.inputs[["Tx_t_max"]], 
    p.Tx                 = l.inputs[["p.Tx"]], 
    Tx                   = l.inputs[["Tx"]],
    v.TX_val             = l.inputs[["v.TX_val"]]
  )
  

  return(a.out)
  
}


#### Module Functions ####
######################################## TX2 (Non-DMT)

f.update_TX2 <- function(v.TX2.lag) {
  tx2 <- v.TX2.lag
  return(tx2)
}


######################################## LTC

f.update_LTC <- function(v.LTC.lag) {
  ltc <- v.LTC.lag
  return(ltc)
}


######################################## TX (DMT)

f.update_TX <- function(v.TX.lag, v.SEV_obs, v.SEV_obs.lag, n.alive, random_cycle, Tx_t_max, p.Tx, Tx, v.TX_val) {
  
  # start with empty vector
  tx <- rep(0, n.alive)
  
  # previous Tx state * previous syndrome==MCI * time shorter than maximum treatment duration
  if(Tx==1) {
    
    initiate_tx <- v.SEV_obs.lag == -9 & (v.SEV_obs == 0 | v.SEV_obs == 1)
    tx[initiate_tx] <- f.qcat(p_rand = random_cycle[initiate_tx], p_cat = p.Tx, values = v.TX_val)
    
    continue_tx <- v.TX.lag == 1 & v.SEV_obs < 2
    tx[continue_tx] <- v.TX.lag[continue_tx] * as.numeric(2 < Tx_t_max)
  }
  
  
  return(tx)
  
}

