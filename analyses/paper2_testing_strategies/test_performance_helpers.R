######################################## GRAM HELPER FUNCTIONS: TEST PERFORMANCE ANALYSES ########################################
require(ggplot2)
require(tidyr)
require(abind)


# Run test performance analyses per cycle
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

# Wrapper function for f.analyze_test_performance, to produce reader-friendly tibble output
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

# Plot point-in-time test performance results for multiple scenarios at once
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


# Post-process raw model output to add summary variables of interest
post_processing_outputs <- function(output_array) {
  
  bha_results <- output_array[,"BHA",]
  pcp_results <- output_array[,"PCP",]
  bha_dx_results <- ifelse(!is.na(pcp_results) & pcp_results >= 0, pcp_results, bha_results)
  
  dx_results <- apply(bha_dx_results, 2, cummax)
  
  syn_status <- output_array[, "SYN", ]
  
  # 1. TPs
  tp_direct <- (dx_results == 1) & (syn_status == 1)
  
  # 2. FPs
  ever_fp <- apply(((dx_results == 1) & (syn_status == 0)), 2, cummax)
  
  # Split FPs into early detection vs actual FP.
  will_be_impaired <- colSums(syn_status == 1, na.rm = TRUE) > 0
  
  real_fp <- t(t(ever_fp) & !will_be_impaired)
  early_pos <- t(t(ever_fp) & will_be_impaired) & (syn_status != 1)
  converted_tp <- t(t(ever_fp) & will_be_impaired) & (syn_status == 1)
  
  # TN
  tn <- (dx_results == 0) & (syn_status < 1)
  
  # FN
  fn <- (dx_results == 0) & (syn_status == 1)
  
  # eligible non-testers
  
  clinical_dx <- output_array[,"DX",]  # DX is carried over in the base model, no need for cummax
  has_hcare <- output_array[,"HCARE",]
  eligible <- clinical_dx == 0 & has_hcare == 1
  
  ever_tested <- apply(bha_results,2,cummax) >= 0
  
  notest_tn <- eligible & !ever_tested & syn_status != 1
  notest_fn <- eligible & !ever_tested & syn_status == 1
  
  # deaths
  deaths <- output_array[,"ALIVE",] == 0
  
  # create data frame
  plot_data <- data.frame(
    age = 50:100,
    tp = rowSums(tp_direct, na.rm = TRUE),
    fp = rowSums(real_fp, na.rm = TRUE),
    early_pos = rowSums(early_pos, na.rm = TRUE),
    converted_tp = rowSums(converted_tp, na.rm = TRUE),
    tn = rowSums(tn, na.rm = TRUE),
    fn = rowSums(fn, na.rm = TRUE),
    notest_tn = rowSums(notest_tn, na.rm = TRUE),
    notest_fn = rowSums(notest_fn, na.rm = TRUE),
    clinical_dx = rowSums(clinical_dx, na.rm = TRUE),
    death = rowSums(deaths))
  
  return(plot_data)
}


# Plot cumulative counts for one or more variables
plot_cumulative_count <- function(output_array, 
                                  variables, # a list of named variables, each with a lower-level list with variable_name and condition_value
                                  y_axis_label = "Count",
                                  cycles = 1:dim(output_array)[1], 
                                  plot_title = "Cumulative Plot",
                                  scenario_name = NULL) {
  
  all_results <- list()
  
  # Loop over each variable condition provided by the user
  for (metric_name in names(variables)) {
    var_info <- variables[[metric_name]]
    variable_name <- var_info$variable_name
    condition_value <- var_info$condition_value
    
    n_individuals <- dim(output_array)[3]
    has_condition <- rep(FALSE, n_individuals)
    cumulative_counts <- integer(length(cycles))
    
    # Calculate cumulative counts for the current variable
    for (i in seq_along(cycles)) {
      t <- cycles[i]
      current_values <- output_array[t, variable_name, ]
      newly_met_condition <- which(current_values %in% condition_value & !has_condition)
      if (length(newly_met_condition) > 0) {
        has_condition[newly_met_condition] <- TRUE
      }
      
      cumulative_counts[i] <- sum(has_condition)
    }
    
    # Store results in a data frame
    all_results[[metric_name]] <- data.frame(
      Age = (50 + cycles) - 1,
      Count = cumulative_counts,
      Metric = metric_name
    )
    
  }
  
  # Combine all data frames into one for plotting
  results_df <- do.call(rbind, all_results)
  
  # Plot with multiple lines
  p <- ggplot(results_df, aes(x = Age, y = Count, color = Metric, group = Metric)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    labs(title = plot_title,
         subtitle = scenario_name,
         x = "Age",
         y = y_axis_label,
         color = "Metric") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = results_df$Age[seq(1, length(unique(results_df$Age)), by = 2)])
  
  return(p)
}


# Plot all the ways to be positive 
plot_test_results <- function(plot_data, 
                              y_axis_label = "Count",
                              ages = 50:100, 
                              plot_title = "Testing Strategy Results",
                              scenario_names = NULL, #must be a named list, where names correspond to scenario value in data
                              y_max = 100000,
                              y_transform = NULL,
                              show_early_pos = TRUE,
                              show_non_testers = TRUE) {
  
  if(!show_early_pos) {
    plot_data <- plot_data %>%
      mutate(fp = fp + early_pos,
             tp = tp + converted_tp) %>%
      select(-early_pos, -converted_tp)
  }
  
  if(!show_non_testers) {
    plot_data <- plot_data %>%
      select(-notest_tn, - notest_fn)
  }
  
  category_map <- tribble(
    ~Result,        ~Status,       ~Test_Value,
    "tp",           "Impaired",    "Positive",
    "fp",           "Healthy",     "Positive",
    "early_pos",    "Healthy",     "Positive",
    "converted_tp", "Impaired",    "Positive",
    "tn",           "Healthy",     "Negative",
    "fn",           "Impaired",    "Negative",
    "notest_tn",    "Healthy",     "No Test",
    "notest_fn",    "Impaired",    "No Test",
    "death",        "Deaths",      NA)
    
    
  plot_data <- plot_data %>%
    select(-clinical_dx) %>%
    filter(age %in% ages) %>%
    mutate(tn_zero_rank = cumsum(notest_tn == 0),
           notest_tn = ifelse((notest_tn == 0) & (tn_zero_rank > 1), 
                              NA, notest_tn)) %>%
    mutate(fn_zero_rank = cumsum(notest_fn == 0),
           notest_fn = ifelse((notest_fn == 0) & (fn_zero_rank > 1), 
                              NA, notest_fn)) %>%
    select(-tn_zero_rank, -fn_zero_rank) %>%
    pivot_longer(cols = c(-age,-scenario),
                 names_to = "Result",
                 values_to = "Result_Value") %>%
    mutate(Result = factor(Result)) %>%
    left_join(category_map, by = "Result") %>%
    rename(Age = age, Scenario = scenario)
  
  no_test_data <- plot_data %>%
    filter(Test_Value == "No Test") %>%
    mutate(no_test_zeros = cumsum(Test_Value)) 
  
  
  # colors for true status (red = impaired, blue = healthy)
  my_colors <- c("Impaired" = "firebrick",
                 "Healthy"  = "deepskyblue",
                 "Deaths"   = "grey")
  
  
  # markers (shape) for test result (21 (filled) for positive, 1 (empty) for negative)
   my_fills <- c(
    "Positive" = "firebrick",      
    "Negative" = "deepskyblue",
    "No Test" = NA)    
  
  my_shapes <- c(
    "Positive" = 21, # Filled circle (21) for Positive and Negative
    "Negative" = 21,
    "No Test" = 4)
  
  
  # linetypes for disease and test concordance (solid doesn't switch between lines, dot-dash does)
  my_linetypes <- c("tp"        = "solid", 
                    "fp"        = "solid", 
                    "early_pos" = "dotted", 
                    "converted_tp" = "dotted",
                    "tn"        = "solid",
                    "fn"        = "solid",
                    "notest_tn" = "solid",
                    "notest_fn" = "solid",
                    "death"     = "solid")
  
  my_linetype_labels <- c(
    "tp"        = "Concordent", 
    "fp"        = "Concordent", 
    "early_pos" = "Discordant", 
    "converted_tp" = "Discordant",
    "tn"        = "Concordent",
    "fn"        = "Concordent",
    "notest_tn" = "Concordent",
    "notest_fn" = "Concordent",
    "death"     = "Concordent"
  )
  
  
  p <- ggplot(plot_data, aes(x = Age, y = Result_Value)) +
    geom_line(aes(group = Result, color = Status), linewidth = 1.5) +
    geom_point(data = subset(plot_data, Test_Value != "No Test"),
               aes(fill = Test_Value,
                   shape = Test_Value),
               color = "white",
               size = 3) +
    geom_point(data = subset(plot_data, Test_Value == "No Test"),
               aes(shape = Test_Value),
               color = "grey10",
               size = 3) +
    facet_wrap(~Scenario, labeller = as_labeller(scenario_names)) + 
    labs(title = plot_title,
         x = "Age",
         y = "Count") +
    ylim(0, y_max) + 
    scale_color_manual(values = my_colors, breaks = names(my_colors)) +    
    scale_shape_manual(values = my_shapes, breaks = names(my_shapes)) +
    scale_fill_manual(values = my_fills, breaks = names(my_fills)) +
    guides(color = guide_legend(title = "Cognitive Status",
                                override.aes = list(size = 3, linewidth = 1.5)),
           shape = guide_legend(title = "Test Result",
                                override.aes = list(fill = my_fills,
                                                    size = 3)),
           fill = "none") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), 
          legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 10)) 
  
  if (!is.null(y_transform)) {
    p <- p + scale_y_continuous(trans = y_transform)
  }
  
  return(p)
  
}


plot_testers <- function(plot_data, 
                         y_axis_label = "Count",
                         ages = 50:100, 
                         plot_title = "Testing Status",
                         scenario_names = NULL, #must be a named list, where names correspond to scenario value in data
                         y_max = 100000) {
  
  
  category_map <- tribble(
    ~Result,        ~Status,       ~Test_Value,
    "all_testers",  NA,            "Tested",
    "notest_tn",    "Healthy",     "Not Tested",
    "notest_fn",    "Impaired",    "Not Tested",
    "death",        "Deaths",      NA)
  
  
  plot_data <- plot_data %>%
    select(-clinical_dx) %>%
    filter(age %in% ages) %>%
    mutate(all_testers = tp + fp + tn + fn + early_pos + converted_tp,
           .keep = "unused") %>%
    mutate(tn_zero_rank = cumsum(notest_tn == 0),
           notest_tn = ifelse((notest_tn == 0) & (tn_zero_rank > 1), 
                              NA, notest_tn)) %>%
    mutate(fn_zero_rank = cumsum(notest_fn == 0),
           notest_fn = ifelse((notest_fn == 0) & (fn_zero_rank > 1), 
                              NA, notest_fn)) %>%
    select(-tn_zero_rank, -fn_zero_rank) %>%
    pivot_longer(cols = c(-age,-scenario),
                 names_to = "Result",
                 values_to = "Result_Value") %>%
    mutate(Result = factor(Result)) %>%
    left_join(category_map, by = "Result") %>%
    rename(Age = age, Scenario = scenario)
  
  # colors for true status (red = impaired, blue = healthy)
  my_colors <- c("Impaired" = "firebrick",
                 "Healthy"  = "deepskyblue",
                 "Deaths"   = "grey")
  
  my_shapes <- c(
    "Tested" = 16, 
    "Not Tested" = 4)
  
  
  p <- ggplot(plot_data, aes(x = Age, y = Result_Value)) +
    geom_line(aes(group = Result, color = Status)) +
    geom_point(aes(shape = Test_Value)) +
    facet_wrap(~fct_rev(Scenario), labeller = as_labeller(scenario_names)) + 
    labs(title = plot_title,
         x = "Age",
         y = "Count") +
    ylim(0, y_max) + 
    scale_color_manual(values = my_colors, breaks = names(my_colors)) +    
    scale_shape_manual(values = my_shapes, breaks = names(my_shapes)) +
    guides(color = guide_legend(title = "Cognitive Status",
                                override.aes = list(size = 3, linewidth = 1.5)),
           shape = guide_legend(title = "Test Status",
                                override.aes = list(size = 3))) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), 
          legend.position = "bottom",
          legend.box = "vertical",
          strip.text = element_text(size = 10)) 
  
  return(p)
    
  
  
  
}