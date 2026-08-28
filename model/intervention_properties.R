######################################## GRAM: INTERVENTION PROPERTIES ########################################
# Treatment and prevention parameters, and the scaffolding the two-arm comparison wrapper
# (f.out_summary) fills in. Kept out of model/setup.R so that setup.R describes natural
# history only.
#
# Sourced by model/helpers/source_all.R, i.e. after model/setup.R has created l.inputs.


## Strategy scaffolding
# Tx is read by f.update_TX(). The remaining five are currently read nowhere in the
# repository; they are the placeholders f.out_summary() was written around.
l.inputs[["Tx"]] <- 0                                      # empty parameter to be filled in as part of the strategies
l.inputs[["strategy"]] <- NA                               # empty parameter to be filled in as part of the strategies
l.inputs[["strategy_strat1"]] <- "control"
l.inputs[["strategy_strat2"]] <- "intervention_dmt"
l.inputs[["Tx_strat1"]] <- 0
l.inputs[["Tx_strat2"]] <- 1


## Treatment effect and eligibility
l.inputs[["rr.Tx_mci"]] <- 0.70
l.inputs[["Tx_t_max"]] <- 3
l.inputs[["p.Tx"]] <- c(0,1) # Probability of DMT ineligible, vs. eligible
l.inputs[["rr.Px_mci"]] <- 1 # Hypothetical -- risk ratio for developing MCI given a prevention intervention (effectiveness of intervention)


## Treatment costs
l.inputs[["c.Tx"]] <- 5000   # cost of DMT
l.inputs[["c.Tx2"]] <- 500   # cost of non-DMT treatment
