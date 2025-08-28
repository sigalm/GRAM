######################################## GRAM HELPER FUNCTIONS: TEST PERFORMANCE ANALYSES ########################################



f.analyze_test_performance <- function(scenario, cycle) {
  
  concordance_vector <- scenario$aggregated_results_totpop$state_concordance[cycle, 2:13]
  
  states <- c("h", "tci", "mci", "dem")
  test_results <- c("neg","pos")
  concordance_table <- as.data.frame(matrix(0, nrow = 4, ncol = 3, 
                                            dimnames = list(states, c(test_results, "NA"))))
  
  for (i in (seq_along(concordance_vector))) {
    split_name <- strsplit(names(concordance_vector)[i], "_")[[1]]
    origin <- split_name[1]
    dest <- split_name[2]
    concordance_table[origin, dest] <- concordance_vector[i]
  }
  
  concordance_table <- concordance_table %>%
    rename(BHApos = pos,
           BHAneg = neg,
           NoTest = `NA`) %>%
    mutate(RowTot =  rowSums(.))
  
  concordance_table_pct <- concordance_table %>%
    mutate(across(-RowTot, ~ .x / RowTot)) %>%
    select(-RowTot) %>%
    {round(.,4)}
  
  
  total_tests <- sum(concordance_table$RowTot) - sum(concordance_table$NoTest) 
  
  specificity <- sum(concordance_table[c("h","tci"),"BHAneg"]) / sum(concordance_table[c("h","tci"),c("BHAneg","BHApos")])
  sensitivity <- c(mci = concordance_table["mci","BHApos"] / sum(concordance_table["mci",c("BHAneg","BHApos")]),
                   dem = concordance_table["dem","BHApos"] / sum(concordance_table["dem",c("BHAneg","BHApos")]))
  ppv <- sum(concordance_table[c("mci","dem"), "BHApos"]) / sum(concordance_table[,"BHApos"])
  npv_test <- sum(concordance_table[c("h","tci"),"BHAneg"]) / sum(concordance_table[,"BHAneg"])
  npv_strat <- sum(concordance_table[c("h","tci"),c("BHAneg","NoTest")]) / sum(concordance_table[,c("BHAneg","NoTest")])
  
  missed_pct <- c(mci_missed = round(sum(concordance_table_pct["mci",c("BHAneg","NoTest")]) * 100, 2),
                  dem_missed = round(sum(concordance_table_pct["dem",c("BHAneg","NoTest")]) * 100, 2))
  
  tci_fu <- round(concordance_table_pct["tci","BHApos"] * 100, 2)
  
  
  return(list(
    total_tests = total_tests,
    concordance_table = concordance_table,
    concordance_table_pct = concordance_table_pct,
    sens = sensitivity, spec = specificity,
    ppv = ppv, npv_test = npv_test, npv_strat = npv_strat,
    missed_pct = missed_pct,
    tci_fu = tci_fu
  ))
}


compute_results <- function(scenario) {
  map_dfr(50:100, function(age) {
    testperf <- f.analyze_test_performance(scenario, age - 50 + 1)
    tibble(
      Scenario = scenario$inputs[["scenario"]][["title"]],
      Age = age,
      N_Tests = testperf$total_tests,
      TP = sum(testperf$concordance_table[c("mci","dem"), "BHApos"]),
      FN = sum(testperf$concordance_table[c("mci","dem"), "BHAneg"]),
      TN = sum(testperf$concordance_table[c("h","tci"), "BHAneg"]),
      FP = sum(testperf$concordance_table[c("h","tci"), "BHApos"]),
      FN_All = sum(testperf$concordance_table[c("mci","dem"), c("BHAneg","NoTest")]),
      TN_All = sum(testperf$concordance_table[c("h","tci"), c("BHAneg","NoTest")]),
      PPV = round(testperf$ppv, 3),
      NPV_Test = round(testperf$npv_test, 3),
      NPV_Strategy = round(testperf$npv_strat, 3)
    )
  })
}


plot_results <- function(scenario_list) {

  all_results <- map_dfr(scenario_list, compute_results)
  
  pv_results <- all_results %>%
    pivot_longer(cols = c(PPV, NPV_Test, NPV_Strategy), names_to = "Metric", values_to = "Value") %>%
    select(Scenario, Age, Metric, Value)
  
  fig_results <- ggplot(data = pv_results, aes(x = Age, y = Value, color = Scenario, linetype = Scenario)) +
    geom_smooth(se = FALSE) + 
    scale_color_manual(values = scenario_colors) +
    scale_linetype_manual(values = scenario_linetypes) +
    xlim(NA, 90) +
    ylim(0,1) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right",
          legend.direction = "vertical")

  fig_results
  
  return(list(results_table = all_results,
              results_figure = fig_results))
}


#' Plot cumulative PPV/NPV results by Age for each scenario
#' Identical to plot_results, but uses cumulative PPV/NPV from cumulative counts
plot_results_cumulative <- function(inputs_list, scenario_list) {
  all_results <- map2_dfr(inputs_list, scenario_list, compute_results)

  # Calculate cumulative sums for TP, FP, TN, FN by Scenario
  cum_results <- all_results %>%
    dplyr::group_by(Scenario) %>%
    dplyr::arrange(Age, .by_group = TRUE) %>%
    dplyr::mutate(
      cum_TP = cumsum(TP),
      cum_FP = cumsum(FP),
      cum_TN = cumsum(TN),
      cum_FN = cumsum(FN),
      cum_TN_All = cumsum(TN_All),
      cum_FN_All = cumsum(FN_All),
      cum_PPV = ifelse(cum_TP + cum_FP > 0, cum_TP / (cum_TP + cum_FP), NA_real_),
      cum_NPV_Test = ifelse(cum_TN + cum_FN > 0, cum_TN / (cum_TN + cum_FN), NA_real_),
      cum_NPV_Strategy = ifelse(cum_TN_All + cum_FN_All > 0, cum_TN_All / (cum_TN_All + cum_FN_All), NA_real_)
    ) %>%
    dplyr::ungroup()

  pv_results_cum <- cum_results %>%
    select(Scenario, Age, cum_PPV, cum_NPV_Test, cum_NPV_Strategy) %>%
    tidyr::pivot_longer(cols = c(cum_PPV, cum_NPV_Test, cum_NPV_Strategy), names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(Metric = dplyr::recode(Metric,
                                         cum_PPV = "PPV",
                                         cum_NPV_Test = "NPV_Test",
                                         cum_NPV_Strategy = "NPV_Strategy"))

  fig_results <- ggplot(data = pv_results_cum, aes(x = Age, y = Value, color = Scenario, linetype = Scenario)) +
    geom_smooth(se = FALSE) + 
    scale_color_manual(values = scenario_colors) +
    scale_linetype_manual(values = scenario_linetypes) +
    xlim(NA, 90) +
    facet_wrap(~ Metric) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "right",
          legend.direction = "vertical")

  fig_results

  return(list(results_table = all_results,
              results_figure = fig_results))
}

