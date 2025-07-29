######################################## GRAM MODULE: MORTALITY ########################################
# This section defines the functions for updating mortality and alive status variables.

f.module_mortality <- function(a.out, t, a.random) {
  
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

  
  return(a.out)
  
}


######################################## ALIVE

# https://adv-r.hadley.nz/subsetting.html vectorized subset from array (i.e. create matrix with array coordinates in each row to be looked up) 

f.update_ALIVE <- function(
    alive.lag,
    v.AGE.lag, v.SYN.lag, v.SEV.lag, 
    m.lifetable, 
    hr.mort_mci, hr.mort_mil, hr.mort_mod, hr.mort_sev, 
    hr.mort_mci_age, hr.mort_mil_age, hr.mort_mod_age, hr.mort_sev_age,
    random_cycle) 
{
  
  # generate life table coordinates for looking up age-specific mortality 
  # (see https://adv-r.hadley.nz/subsetting.html paragraph 4.2.3 subsetting > selecting multiple elements > subsetting)
  lifetable_lookup_coordinates <- matrix(data = round(v.AGE.lag,0), ncol = 1) - 50 + 1
  
  # determine relative mortality risk related to syndrome and severity
  healthy <- v.SYN.lag<1
  mci <- v.SYN.lag==1 & v.SEV.lag==0
  mil <- v.SYN.lag==1 & v.SEV.lag==1
  mod <- v.SYN.lag==1 & v.SEV.lag==2
  sev <- v.SYN.lag==1 & v.SEV.lag==3
  
  v.risk_ratio <- rep(NA, length(v.AGE.lag))
  v.risk_ratio[healthy] <- 1
  v.risk_ratio[mci] <- hr.mort_mci * ((v.AGE.lag[mci]<70) + hr.mort_mci_age[2] * (v.AGE.lag[mci]>=70 & v.AGE.lag[mci]<80) + hr.mort_mci_age[3] * (v.AGE.lag[mci]>=80))
  v.risk_ratio[mil] <- hr.mort_mil * ((v.AGE.lag[mil]<70) + hr.mort_mil_age[2] * (v.AGE.lag[mil]>=70 & v.AGE.lag[mil]<80) + hr.mort_mil_age[3] * (v.AGE.lag[mil]>=80))
  v.risk_ratio[mod] <- hr.mort_mod * ((v.AGE.lag[mod]<70) + hr.mort_mod_age[2] * (v.AGE.lag[mod]>=70 & v.AGE.lag[mod]<80) + hr.mort_mod_age[3] * (v.AGE.lag[mod]>=80))
  v.risk_ratio[sev] <- hr.mort_sev * ((v.AGE.lag[sev]<70) + hr.mort_sev_age[2] * (v.AGE.lag[sev]>=70 & v.AGE.lag[sev]<80) + hr.mort_sev_age[3] * (v.AGE.lag[sev]>=80))
  
  
  # calculate death probability using life table, relative risk and adjustment for cycle length
  prob_death <- pmin(f.adjustprobability(
    p = m.lifetable[lifetable_lookup_coordinates], 
    t_new = 1, 
    t_old = 1,
    RR = v.risk_ratio),1)
  
  # compare probability to random value
  alive <- as.numeric(!prob_death > random_cycle)
  
  # check for errors
  if (is.na(sum(alive))) {
    print(table(alive, useNA = "always"))
    print(table(alive, v.SYN.lag, useNA = "always"))
    print(table(alive, v.SEV.lag, useNA = "always"))
    stop("NAs produced in 'ALIVE'")
  }
  # return
  return(alive)
}

