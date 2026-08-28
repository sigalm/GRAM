######################################## GRAM: COGNITIVE TEST PROPERTIES ########################################
# Properties of the cognitive tests themselves: what an instrument does when it is
# administered, and what one administration costs. These are measurement characteristics
# of the test, not choices about a testing programme -- WHO is tested, HOW OFTEN, and
# WHAT FOLLOWS a positive result are defined per scenario in
# analyses/*/bha_scenarios/*_config.R and belong nowhere else.
#
# Sourced by model/helpers/source_all.R, i.e. after model/setup.R has created l.inputs.


## Test performance
# Source: Elena Tsoy (both CS and GS at the -1.5z cutoff)
#
# The sensitivity vector is indexed by the branch order in f.update_BHA(), NOT by severity:
#   [1] TCI (SYN == 0.5)
#   [2] non-progressive memory loss (MEMLOSS == 1)
#   [3] MCI (SEV == 0)
#   [4] dementia (SEV >= 1 -- all three dementia severities share this one value)
#
# Exposed as plain objects rather than l.inputs entries so that a scenario config can
# reference them directly. Configs are sourced into an environment whose parent chain
# reaches the global environment, so an l.inputs lookup inside a config resolves to the
# global l.inputs and silently ignores whatever inputs object was passed to
# load_scenario(). Referring to BHA_GS avoids that trap.
BHA_CS <- list(sens = c(0.44, 0.48, 0.66, 0.96), spec = 0.93)
BHA_GS <- list(sens = c(0.56, 0.61, 0.84, 0.98), spec = 0.92)


## Rater error
l.inputs[["r.CDR_sd3"]] <- 0   # rater error (inter-rater reliability, will be less reliable in MCI, better in dem)


## Unit costs of the testing pathway
l.inputs[["c.bha"]] <- 200       # cost of administering the BHA
l.inputs[["c.bhapos"]] <- 2000   # cost of follow up with patient with positive BHA
l.inputs[["c.np"]] <- 1000       # cost of neuropsych assessment
