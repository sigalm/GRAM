######################################## GRAM HELPER FUNCTIONS: TEST PERFORMANCE ANALYSES ########################################



f.analyze_test_performance <- function(scenario, cycle) {
  
  concordance_vector <- scenario$aggregated_results_totpop$state_concordance[cycle, 1:12] * ncol(scenario$output[cycle,,])
  n_alive <- sum(concordance_vector)
  
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
  npv <- sum(concordance_table[c("h","tci"),"BHAneg"]) / sum(concordance_table[,"BHAneg"])
  
  missed_pct <- c(mci_missed = round(sum(concordance_table_pct["mci",c("BHAneg","NoTest")]) * 100, 2),
                  dem_missed = round(sum(concordance_table_pct["dem",c("BHAneg","NoTest")]) * 100, 2))
  
  tci_fu <- round(concordance_table_pct["tci","BHApos"] * 100, 2)
  
  
  return(list(
    total_tests = total_tests,
    concordance_table = concordance_table,
    concordance_table_pct = concordance_table_pct,
    sens = sensitivity, spec = specificity,
    ppv = ppv, npv = npv,
    missed_pct = missed_pct,
    tci_fu = tci_fu
  ))
}
