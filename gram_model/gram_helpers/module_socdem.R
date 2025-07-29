######################################## GRAM MODULE: SOCIODEMPGRAPHIC ATTRIBUTES ########################################
# This section defines the functions for updating sociodemographic attributes.

#### Module Wrapper ####
f.module_socdem <- function(a.out, t, a.random, alive) {
  
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
  
  # HCARE
  a.out[t,"HCARE",alive] <- f.update_HCARE(
    v.HCARE.lag  = a.out[t-1,"HCARE",alive],
    v.AGE        = a.out[t,"AGE",alive],
    random_cycle = a.random[t,"HCARE",alive]
  )    
  
  return(a.out)
  
}


#### Module Functions ####
######################################## AGE

f.update_AGE <- function(v.AGE.lag) {
  age <- v.AGE.lag + 1 # 1 represents fixed cycle time of 1 year
  return(age)
}

######################################## SEX

f.update_SEX <- function(v.SEX.lag) {
  sex <- v.SEX.lag
  return(sex)
}

######################################## RACEETH

f.update_RACEETH <- function(v.RACEETH.lag) {
  raceeth <- v.RACEETH.lag
  return(raceeth)
}

######################################## INCOME

f.update_INCOME <- function(v.INCOME.lag) {
  income <- v.INCOME.lag
  return(income)
}

######################################## EDU

f.update_EDU <- function(v.EDU.lag) {
  edu <- v.EDU.lag
  return(edu)
}

######################################## HCARE

f.update_HCARE <- function(v.HCARE.lag, v.AGE, random_cycle) {
  hcare <- v.HCARE.lag
  
  hcare[v.AGE == 65 & v.HCARE.lag == 0] <- f.qcat(p_rand = random_cycle[v.AGE == 65 & v.HCARE.lag == 0],
                                                  p_cat = c(0.15, 0.85), values = c(0,1))
  
  return(hcare)
}