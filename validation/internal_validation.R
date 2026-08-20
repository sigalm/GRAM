# ****************************************************************
# ======= GRAM INTERNAL VALIDATION: WHICAP =======
# ****************************************************************

# Run GRAM using the WHICAP cohort as the starting population. This is the cohort in which the MCI hazard function was estimated.
# Goal is to replicate MCI incidence from that same study: Angevaare et al 2022.

# The general model setup is sourced rather than copied, so this script only states what differs
# from it. No microdata is supplied, so the cohort is built from defaults and the parameters below.

source("model/setup.R") # note: this clears the environment, so keep it first
source("model/helpers/source_all.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")

library(tableone)
library(easystats)


######################################## COHORT OVERRIDES ########################################

## Model settings
l.inputs1 <- l.inputs
l.inputs1[["n.ind"]] <- 2000                         # number of individuals to simulate
l.inputs1[["n.cycle"]] <- 12                          # number of cycles to simulate (2 more than Framingham because last cycle of model is flat and 1st cycle is baseline)
l.inputs1[["seed_stochastic"]] <- 20250624

## Demographics of the WHICAP cohort
l.inputs1[["AGE_start_mean"]] <- 73                  
l.inputs1[["AGE_start_sd"]] <- 6
l.inputs1[["p.SEX_start_male"]] <- 1-0.669
l.inputs1[["p.SEX_start_female"]] <- 0.669
l.inputs1[["p.APOE4_start"]] <- c(1-0.263, 0.263)        # p for non-carrier and carrier, respectively

# << need to change to accept distn as continuous var >> 
l.inputs1[["EDU_start_mean"]]   <- 10.9
l.inputs1[["EDU_start_sd"]] <- 5

l.inputs1[["v.RACEETH_val"]] <- c(0, 1, 2)    # 0 = White, 1 = Black, 2 = Hispanic
l.inputs1[["p.RACEETH_start"]] <- c((1-(0.293+0.464)), 0.293, 0.464) 

l.inputs1[["MEDBUR_start_mean"]] <- 2.40
l.inputs1[["MEDBUR_start_sd"]] <- 1.6

l.inputs1[["v.INCOME_val"]]  <- c(0, 1, 2)    # 0 = low (<$9000/y), 1 = medium ($9000-$36000/y), 2 = high (>$36000/y)
l.inputs1[["p.INCOME_start"]]  <- c(0.261, 0.558, 0.181) 

l.inputs1[["p.SYN_start"]] <- c(1, 0, 0)    # p for SYN == 0 (normal), 0.5 (TCI) and 1 (impaired)
l.inputs1[["p.MEMLOSS_start"]] <- c(0.91, 0.09)

## Mortality
l.inputs1[["hr.mort_mci_age"]] <- c(1, 1, 1)
l.inputs1[["hr.mort_mod_age"]] <- l.inputs[["hr.mort_sev_age"]] <- c(1, 1, 1)
l.inputs1[["m.lifetable"]][nrow(l.inputs[["m.lifetable"]]), ] <- 1   # certain death at the oldest tabulated age

## Grab calibrated parameters
l.inputs1[["r.CDRslow_mean"]] <-  
  (seq(0, 1, length.out = 51)^1.9500) * (1.8997 * l.inputs[["r.CDRslow_mean"]])
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 1.6408


######################################## RUN ########################################

sim_int_val <- f.run(l.inputs1, microdata = NULL, printLevel = 0)

# << move this to the epi_helpers later >>
calculate_mci_incidence <- function(a.out) {
  new_cases <- 0
  person_years <- 0
  for (t in 4:(dim(a.out)[1]-1)) {  # 1 less than cycle count because last two cycles are flat
    # Select individuals in the age group who were healthy at the start of the cycle
    at_risk <- which(a.out[t - 1, "SYN", ] != 1)
    
    # Calculate person-years at risk for this time point
    person_years <- person_years + length(at_risk)
    
    # Count new MCI cases (healthy in previous cycle and MCI in current cycle)
    new_cases <- new_cases + sum(a.out[t, "SEV", at_risk] == 0, na.rm = TRUE)
  }
  incidence_rate <- if (person_years > 0) new_cases / person_years else NA
  
  ci <- poisson.test(new_cases, person_years)$conf.int
  
  return(
    list(incident_cases = new_cases,
         person_years = person_years,
         incidence_rate = incidence_rate,
         i_lower       = ci[1],
         ci_upper       = ci[2]
         ))
}

int_val_incidence <- calculate_mci_incidence(sim_int_val)

# Angevaare et al. reported ~3 fewer years of follow-up, on average, for participants who did
# not develop incident MCI compared to those who did. This reflects non-mortality loss to
# follow-up (missed visits, dropout, staggered study end) that GRAM has no mechanism to
# reproduce -- every simulated individual is followed until death or the end of the window.
# Individuals who died during follow-up already have correctly truncated person-time from
# GRAM's own mortality module, so the correction is restricted to non-cases who survived the
# full 8-cycle follow-up window; subtracting 3 more years from someone whose time was already
# cut short by simulated death would double-count censoring.
last_cycle_followup <- dim(sim_int_val)[1] - 1  # cycle 11: last of the 8 follow-up cycles
ever_mci <- apply(sim_int_val[4:last_cycle_followup, "SEV", ] == 0, 2, any, na.rm = TRUE)
cases <- which(ever_mci)

alive_baseline <- which(sim_int_val[3, "ALIVE", ] == 1)  # cycle 3: baseline visit
noncases_survived_followup <- intersect(alive_baseline, which(!ever_mci))
noncases_survived_followup <- noncases_survived_followup[
  sim_int_val[last_cycle_followup, "ALIVE", noncases_survived_followup] == 1
]

person_years_adj <- int_val_incidence$person_years - 3 * length(noncases_survived_followup)
incidence_rate_adj <- int_val_incidence$incident_cases / person_years_adj
incidence_rate_adj * 1000 # per 1,000 person-years
poisson.test(int_val_incidence$incident_cases, person_years_adj)$conf.int * 1000

# 50.13 (95% CI 45.55 - 55.05) in GRAM, compared to Angevaare et al.'s 56 (95% CI: 52 - 60), per 1,000 person years.


baseline_incident_mci <- as.data.frame(t(sim_int_val[3, , cases]))
mean(baseline_incident_mci$AGE)
prop.table(table(baseline_incident_mci$RACEETH))
mean(baseline_incident_mci$MEDBUR)
prop.table(table(baseline_incident_mci$INCOME))
sum(sim_int_val[3,"ALIVE",]==1) # alive at baseline (cycle 3)

