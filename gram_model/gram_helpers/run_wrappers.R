######################################## GRAM MAIN FUNCTIONS: WRAPPER FUNCTIONS  ########################################

######################################### BASIC WRAPPER 
# Creates 3D output, aggregates results, and charts disease progression.

f.wrap_run <- function(l.inputs, microdata = NULL, printLevel = 0) {
  
  output <- f.run(l.inputs = l.inputs, microdata = microdata, printLevel = printLevel)
  
  # Aggregate results for full cohort
  aggregated_results_totpop <- f.out_aggregate(a.out = output, l.inputs = l.inputs)
  fig.progression <- f.make_figures(l.out = aggregated_results_totpop, l.inputs = l.inputs)
  
  # # Aggregate results for those who develop MCI at any point
  # impaired <- apply(output[,"SYN",], 2, function(x) any(x == 1, na.rm = TRUE))
  # output_impaired <- output[,,impaired]
  # aggregated_results_impaired <- f.out_aggregate(a.out = output_impaired, l.inputs = l.inputs)
  # figures_impaired <- f.figures(l.out = aggregated_results_impaired, l.inputs = l.inputs)
  # 
  # # Aggregate results for only DMT eligibles
  # treated <- apply(output[,"TX",], 2, function(x) any(x == 1, na.rm = TRUE))
  # output_treated <- output[,,treated]
  # aggregated_results_treated <- f.out_aggregate(a.out = output_treated, l.inputs = l.inputs)
  # figures_treated <- f.figures(l.out = aggregated_results_treated, l.inputs = l.inputs)
  # 
  return(list(
    inputs = l.inputs,
    output = output,
    aggregated_results_totpop = aggregated_results_totpop,
    fig.progression = fig.progression
    # ,
    # aggregated_results_impaired = aggregated_results_impaired,
    # figures_impaired = figures_impaired,
    # aggregated_results_treated = aggregated_results_treated,
    # figures_treated = figures_treated
  )
  )
}

######################################## RUN STRATEGIES
# runs 2 strategies and calculate incremental outcomes

f.out_summary <- function(l.inputs1, l.inputs2, printLevel = 0) {
  
  # strategy 1
  l.inputs_strat1 <- l.inputs1
  # l.inputs_strat1[["strategy"]] <- l.inputs[["strategy_strat1"]]
  # l.inputs_strat1[["Tx"]] <- l.inputs[["Tx_strat1"]]
  a.out_strat1 <- f.run(l.inputs = l.inputs_strat1, printLevel = printLevel)
  a.out_qc_strat1 <- f.qaly_cost(a.out = a.out_strat1, l.inputs = l.inputs_strat1)
  out_strat1 <- f.out_aggregate(a.out = a.out_qc_strat1, l.inputs = l.inputs_strat1)
  fig_strat1 <- f.make_figures(l.out = out_strat1, l.inputs = l.inputs1)
  
  
  # strategy 2
  l.inputs_strat2 <- l.inputs2
  # l.inputs_strat2[["strategy"]] <- l.inputs[["strategy_strat2"]]
  # l.inputs_strat2[["Tx"]] <- l.inputs[["Tx_strat2"]]
  a.out_strat2 <- f.run(l.inputs = l.inputs_strat2, printLevel = printLevel)
  a.out_qc_strat2 <- f.qaly_cost(a.out = a.out_strat2, l.inputs = l.inputs_strat2)
  out_strat2 <- f.out_aggregate(a.out = a.out_qc_strat2, l.inputs = l.inputs_strat2)
  fig_strat2 <- f.make_figures(l.out = out_strat2, l.inputs = l.inputs2)
  
  # summary outcomes
  m.out <- matrix(
    data = NA,
    nrow = 5,
    ncol = ncol(out_strat1[["table_sum"]]),
    dimnames = list( c("strategy 1 (cau)","strategy 2 (dmt)","strategy 3 (symptomatic)","incr. strategy 2-1","incr. strategy 3-1"), colnames(out_strat1[["table_sum"]]) )
  )
  m.out[1,] <- out_strat1[["table_sum"]]
  m.out[2,] <- out_strat2[["table_sum"]]
  
  m.out[4,] <- m.out[2,] - m.out[1,]
  m.out[5,] <- m.out[3,] - m.out[1,]
  m.out <- cbind(m.out, ICER = NA, ICER.dis = NA)
  m.out[4,"ICER"] <- (m.out[2,"COST_tot.sum"    ] - m.out[1,"COST_tot.sum"    ]) / (m.out[2,"QALY.sum"    ] - m.out[1,"QALY.sum"    ])
  m.out[4,"ICER.dis"] <- (m.out[2,"COST_tot.dis.sum"] - m.out[1,"COST_tot.dis.sum"]) / (m.out[2,"QALY.dis.sum"] - m.out[1,"QALY.dis.sum"])
  
  return(list(
    l.inputs1 = l.inputs1,
    l.inputs2 = l.inputs2,
    a.out_strat1 = a.out_strat1, a.out_strat2 = a.out_strat2,
    a.out_qc_strat1 = a.out_qc_strat1, a.out_qc_strat2 = a.out_qc_strat2,
    out_strat1 = out_strat1, out_strat2 = out_strat2,
    fig_strat1 = fig_strat1, fig_strat2 = fig_strat2,
    m.out = m.out
  ))
  
}




