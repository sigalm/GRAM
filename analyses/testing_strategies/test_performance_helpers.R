######################################## GRAM HELPER FUNCTIONS: TEST PERFORMANCE ANALYSES ########################################
require(ggplot2)
require(tidyr)
require(dplyr)
require(scales)
require(patchwork)
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
  
  # Everyone holding a standing positive verdict falls into exactly one of these
  # two totals, on their CURRENT status: a positive in someone not yet impaired
  # is a false positive, and becomes a true positive the moment they convert.
  pos      <- dx_results == 1
  tp_total <- pos & (syn_status == 1)
  fp_total <- pos & (syn_status <  1)
  
  # Early detection is a DISPLAY split only, for show_early_pos = TRUE. The bands
  # are carved OUT of the totals above rather than added to them, so they always
  # sum back exactly. Two things went wrong when they were computed independently:
  # converted_tp was a strict subset of tp_direct (ever_fp implies dx_results == 1,
  # which cummax then carries forward), so tp + converted_tp double-counted; and
  # ever_fp keys on syn_status == 0 rather than < 1, so a first positive during
  # TCI landed in neither FP band and vanished from fp + early_pos.
  #
  # ever_early keeps the syn_status == 0 test on purpose: "early" means flagged
  # before any sign at all, not merely before conversion. TCI-first positives are
  # therefore ordinary false positives, and now land in real_fp instead of nowhere.
  ever_early <- apply(((dx_results == 1) & (syn_status == 0)), 2, cummax)
  will_be_impaired <- colSums(syn_status == 1, na.rm = TRUE) > 0
  early_flag <- t(t(ever_early) & will_be_impaired)
  
  converted_tp <- tp_total &  early_flag    # caught early, impaired now
  tp_direct    <- tp_total & !converted_tp  # flagged when already impaired
  early_pos    <- fp_total &  early_flag    # caught early, not yet impaired
  real_fp      <- fp_total & !early_pos     # positive, will not convert
  
  # TN
  tn <- (dx_results == 0) & (syn_status < 1)
  
  # FN
  fn <- (dx_results == 0) & (syn_status == 1)
  
  # eligible non-testers
  
  clinical_dx <- output_array[,"DX",]  # DX is carried over in the base model, no need for cummax
  has_hcare <- output_array[,"HCARE",]
  
  # Eligibility uses LAGGED DX, matching f.update_BHA: within a cycle the modules run
  # BHA before DX, so the testing decision is made against last cycle's diagnosis. Using
  # the current cycle would drop people from the denominator in the same year they were
  # still eligible to be tested. HCARE is current, which also matches the gate.
  dx_lag <- rbind(rep(NA, ncol(clinical_dx)), clinical_dx[-nrow(clinical_dx), , drop = FALSE])
  eligible <- dx_lag == 0 & has_hcare == 1
  
  ever_tested <- apply(bha_results,2,cummax) >= 0
  
  notest_tn <- eligible & !ever_tested & syn_status != 1
  notest_fn <- eligible & !ever_tested & syn_status == 1
  
  # Alive, never tested, and NOT eligible: no provider, or already carrying a
  # clinical diagnosis. Kept apart from notest_tn / notest_fn on purpose -- those
  # two measure a coverage gap the programme could still close, and someone
  # ineligible is not that. Together with the four tested cells, the two notest_
  # cells and the deaths below, this partitions the starting cohort exactly.
  # NA in cycle 1, where DX has no lag, so the partition only closes from cycle 2.
  not_eligible <- (output_array[,"ALIVE",] == 1) & !ever_tested & !eligible
  
  # deaths, and the survivors they are the complement of
  deaths <- output_array[,"ALIVE",] == 0
  alive  <- output_array[,"ALIVE",] == 1
  
  # create data frame
  plot_data <- data.frame(
    age = 50:(50 + dim(output_array)[1] - 1),
    tp = rowSums(tp_direct, na.rm = TRUE),
    fp = rowSums(real_fp, na.rm = TRUE),
    early_pos = rowSums(early_pos, na.rm = TRUE),
    converted_tp = rowSums(converted_tp, na.rm = TRUE),
    tn = rowSums(tn, na.rm = TRUE),
    fn = rowSums(fn, na.rm = TRUE),
    notest_tn = rowSums(notest_tn, na.rm = TRUE),
    notest_fn = rowSums(notest_fn, na.rm = TRUE),
    clinical_dx = rowSums(clinical_dx, na.rm = TRUE),
    not_eligible = rowSums(not_eligible, na.rm = TRUE),
    alive = rowSums(alive, na.rm = TRUE),
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
    geom_line(aes(group = Result, color = Status)) +
    geom_point(data = subset(plot_data, Test_Value != "No Test"),
               aes(fill = Test_Value,
                   shape = Test_Value),
               color = "white") +
    geom_point(data = subset(plot_data, Test_Value == "No Test"),
               aes(shape = Test_Value),
               color = "grey10") +
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

######################################## COUNTS + COHORT-SHARE FIGURE ########################################
#
# plot_counts_and_shares() draws the two-row main figure: counts per cell above,
# and below, the same cells as a share of the WHOLE starting cohort.
#
# Encoding notes, since this deliberately departs from plot_test_results():
#   - One colour per OUTCOME, not one for status and another for test result.
#     The old scheme put true status in the line colour and test result in the
#     point colour, so a false positive (healthy line, positive points) and a
#     false negative (impaired line, negative points) drew as the same red-and-
#     blue object and whichever came last won. Hue still carries status (warm =
#     impaired, cool = healthy) and shade carries whether the test was right.
#   - Square-root y-axis on the counts row. Reactive is an order of magnitude
#     smaller than Inclusive; on a shared linear axis it is flat against zero.
#   - Mortality on the counts row is carried as SURVIVORS by default, not
#     cumulative deaths. Same information, opposite direction of travel: the
#     survivor line sits above the test series rather than rising up through them.
#   - Series are labelled at the line end in every panel, so neither row needs a
#     legend competing for space.
#
# The bottom row's denominator is the STARTING COHORT, not the tested. Closing
# the stack to 100% keeps growth in reach visible as area rather than dividing it
# out, and leaves the untested and deceased remainder on the page. At any age
# every member of the cohort is in exactly one band:
#
#   dead                                             -> death
#   alive, ever tested                               -> TP / FN / TN / FP
#   alive, never tested, ELIGIBLE                    -> notest_tn / notest_fn
#   alive, never tested, NOT eligible                -> not_eligible
#
# "Not eligible" stays its own band rather than folding into "not tested",
# because notest_tn / notest_fn are deliberately eligibility-gated: they measure
# a coverage gap the programme could still close, and someone ineligible is not
# that. post_processing_outputs() supplies all four groups; they are verified to
# sum to the starting cohort exactly.
#
# All the series are cummax stocks, so someone tested at 65 is still counted at
# 80. That is intentional: the claim is cumulative burden ("23% of the cohort
# carries a false positive by 80"), not per-test yield. Per-test yield would be a
# flow built from the raw array instead.

# Full names, spelled out. "FP" is jargon the reader decodes on every glance.
outcome_long <- c(Deaths = "Deaths", Alive = "Remaining alive",
                  TP = "True positive",  FN = "False negative",
                  TN = "True negative",  FP = "False positive",
                  `No test, impaired` = "Not tested, impaired",
                  `No test, healthy`  = "Not tested, healthy",
                  `No test`           = "Not tested",
                  `Not eligible`      = "Not eligible",
                  Deceased            = "Deceased")

# Deaths take the light grey and sit behind everything, so "not tested, healthy"
# moves off grey onto a muted teal -- cool like the other healthy series, but
# desaturated like its brown "not tested, impaired" partner.
pal_outcome <- c(Deaths = "#C2C7CA", Alive = "#8A9BA3",
                 TP = "#9E2A2B", FN = "#E8A33D",
                 TN = "#1B4965", FP = "#5FA8D3",
                 `No test, impaired` = "#8C6D46", `No test, healthy` = "#6F9A94",
                 `No test`           = "#8FA8A2",
                 `Not eligible`      = "#B5AEA4",
                 Deceased            = "#C3C9CC")

# Fills stay pale so the two largest bands do not dominate, but pale fill makes
# illegible label text on white, so labels take a darker shade of the same hue.
# Assign in place: c(pal_outcome, Deceased = ...) would APPEND a second entry
# under the same name and a [[ ]] lookup returns the FIRST match, so the override
# would be silently ignored and the pale fill colour used for the text.
pal_label <- pal_outcome
pal_label[["Not eligible"]] <- "#5F5A50"
pal_label[["Deceased"]]     <- "#4A5459"
pal_label[["Deaths"]]       <- "#4A5459"

outcome_levels <- c("Deaths", "Alive", "TP", "FN", "TN", "FP",
                    "No test, impaired", "No test, healthy")


# Push colliding labels apart by splitting the difference between the two
# involved, relaxing until nothing moves. A single bottom-up pass instead anchors
# the lowest label and shunts the whole stack upward, leaving labels floating
# above the line ends they belong to. `trans` works in the transformed space so
# this behaves on a square-root axis.
f.spread_labels <- function(y, gap, trans = c("identity", "sqrt")) {
  trans <- match.arg(trans)
  f  <- if (trans == "sqrt") sqrt else identity
  fi <- if (trans == "sqrt") function(x) x^2 else identity
  ord <- order(y)
  z <- f(y[ord])
  for (iter in seq_len(80)) {
    moved <- FALSE
    for (i in seq_along(z)[-1]) {
      deficit <- gap - (z[i] - z[i - 1])
      if (deficit > 1e-9) {
        z[i - 1] <- z[i - 1] - deficit / 2
        z[i]     <- z[i]     + deficit / 2
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  out <- numeric(length(y))
  out[ord] <- fi(pmax(z, 0))
  out
}


plot_counts_and_shares <- function(plot_data,
                                   scenario_names,          # named vector, key -> label; its ORDER sets panel order
                                   strategy_stats = NULL,   # scenario / n_people / n_tests / per_person; NULL = name-only strips
                                   ages = 65:80,
                                   top_title    = "Cumulative strategy-level outcomes among eligible and alive",
                                   top_subtitle = NULL,
                                   bottom_title    = "Share of the cohort in each outcome",
                                   bottom_subtitle = "Stack totals 100% of the starting cohort; the tested block at the bottom is programme reach",
                                   bottom_caption  = "* Not eligible: no healthcare provider, or already has a known cognitive impairment",
                                   y_breaks   = c(0, 1000, 5000, 15000, 30000, 60000),
                                   mortality = c("alive", "deaths", "none"),
                                   split_untested = FALSE,  # bottom row: one "Not tested" band, or cut healthy / impaired
                                   base_size  = 15,     # everything else scales off this
                                   label_size = 3.15,   # line-end and band labels
                                   gap_frac   = 0.062,  # minimum label separation, as a share of the axis
                                   right_pad  = NULL,   # gutter for the labels; defaults to fit label_size
                                   heights    = c(1.15, 1)) {

  keys <- names(scenario_names)

  # Counts row can carry the cohort's mortality either way round: the survivors
  # (the default) or the cumulative deaths they are the complement of. Both come
  # straight from post_processing_outputs(), so no cohort size is needed here.
  mortality <- match.arg(mortality)

  # The gutter has to grow with the type, or bigger labels run off the panel.
  # 0.146 per point of label_size is what fits the longest label at the default.
  if (is.null(right_pad)) right_pad <- 0.146 * label_size

  scen_levels <- unname(scenario_names)

  # Every band, wide, for both rows. The counts row drops the cohort remainder
  # bands; the share row keeps them, which is what closes its stack to 100%.
  wide <- plot_data %>%
    filter(age %in% ages, scenario %in% keys) %>%
    transmute(Age      = age,
              Scenario = factor(scenario_names[scenario], levels = scen_levels),
              TP = tp + converted_tp,
              FP = fp + early_pos,
              TN = tn,
              FN = fn,
              `No test, healthy`  = notest_tn,
              `No test, impaired` = notest_fn,
              `Not eligible`      = not_eligible,
              Deceased            = death,
              Deaths              = death,
              Alive               = alive)

  d <- wide %>%
    select(-`Not eligible`, -Deceased,
           -all_of(setdiff(c("Deaths", "Alive"),
                           switch(mortality, deaths = "Deaths", alive = "Alive", none = character(0))))) %>%
    pivot_longer(-c(Age, Scenario), names_to = "Outcome", values_to = "Count") %>%
    mutate(Outcome = factor(Outcome, levels = outcome_levels),
           Not_tested = grepl("^No test", Outcome))

  cohort_series <- c("Deaths", "Alive")

  # Shared so the two panel grids align and 65-80 sits at the same horizontal
  # position in each row.
  age_breaks <- pretty(range(ages), n = 4)
  age_breaks <- age_breaks[age_breaks >= min(ages) & age_breaks <= max(ages)]
  x_shared <- scale_x_continuous(breaks = age_breaks,
                                 expand = expansion(mult = c(0.03, right_pad)))

  # Self-contained: this figure is a patchwork, so do NOT add theme_paper2 to it
  # afterwards. `&` would push theme_paper2's 20pt bold strip.text onto the
  # bottom row's stats line, which is a long string that has to stay small.
  # Size it with base_size instead.
  base <- theme_minimal(base_size = base_size) +
    theme(panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_line(colour = "grey92"),
          panel.spacing.x    = unit(1.9, "lines"),
          axis.title         = element_text(face = "bold"),
          strip.text         = element_text(size = base_size, face = "bold", hjust = 0),
          plot.title         = element_text(face = "bold", size = base_size * 1.15),
          plot.subtitle      = element_text(colour = "grey35", size = base_size * 0.8),
          legend.position    = "none")

  ## Top row: counts, every series labelled at its line end in every panel
  span <- sqrt(max(d$Count)) - sqrt(min(d$Count))
  line_labels <- d %>%
    filter(Age == max(Age)) %>%
    group_by(Scenario) %>%
    mutate(y   = f.spread_labels(Count, gap_frac * span, "sqrt"),
           lab = outcome_long[as.character(Outcome)]) %>%
    ungroup()

  p_top <- ggplot(d, aes(Age, Count, colour = Outcome, group = Outcome)) +
    # drawn first, so the test series sit on top of it rather than under it
    geom_line(data = filter(d, Outcome %in% cohort_series), linewidth = 1.4) +
    geom_line(data = filter(d, !Outcome %in% cohort_series),
              aes(linetype = Not_tested), linewidth = 1.05) +
    geom_text(data = line_labels, aes(x = Age, y = y, label = lab),
              hjust = 0, nudge_x = 0.3, size = label_size, fontface = "bold",
              inherit.aes = FALSE, colour = pal_label[as.character(line_labels$Outcome)]) +
    facet_wrap(~Scenario) +
    scale_colour_manual(values = pal_outcome) +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "21")) +
    x_shared +
    scale_y_sqrt(labels = comma, breaks = y_breaks) +
    labs(title = top_title, subtitle = top_subtitle, y = "Number of people\n(square-root scale)") +
    base

  ## Bottom row: share of the whole cohort, labelled at the right edge.
  # position_stack puts the FIRST factor level on TOP, so the level order below
  # is the reverse of the bottom-up reading order.
  w_share <- wide %>% select(-Deaths, -Alive)
  if (!split_untested) {
    w_share <- w_share %>%
      mutate(`No test` = `No test, healthy` + `No test, impaired`, .keep = "unused")
  }
  untested_lv <- if (split_untested) c("No test, impaired", "No test, healthy") else "No test"
  band_levels <- c("Deceased", "Not eligible", untested_lv, "FP", "TN", "FN", "TP")

  d_share <- w_share %>%
    pivot_longer(-c(Age, Scenario), names_to = "Band", values_to = "Count") %>%
    mutate(Band = factor(Band, levels = band_levels)) %>%
    group_by(Scenario, Age) %>%
    mutate(Share = Count / sum(Count)) %>%   # sums to exactly 1: verified partition
    ungroup()

  band_labels <- d_share %>%
    filter(Age == max(Age)) %>%
    group_by(Scenario) %>%
    arrange(desc(Band), .by_group = TRUE) %>%
    mutate(ymid = cumsum(Share) - Share / 2,
           y    = f.spread_labels(ymid, 0.075),
           lab  = outcome_long[as.character(Band)],
           # asterisk points at the caption; "not eligible" is the one band whose
           # membership rule is not obvious from its name
           lab  = ifelse(Band == "Not eligible", paste0(lab, "*"), lab)) %>%
    ungroup()

  strip_fn <- if (is.null(strategy_stats)) {
    identity
  } else {
    st <- strategy_stats
    st$label <- unname(scenario_names[st$scenario])
    as_labeller(setNames(
      sprintf("N ever tested = %s\nN total tests done = %s\nTests per person tested = %.1f",
              comma(st$n_people), comma(st$n_tests), st$per_person),
      st$label))
  }

  p_bottom <- ggplot(d_share, aes(Age, Share, fill = Band)) +
    geom_area(colour = "white", linewidth = 0.25) +
    geom_text(data = band_labels, aes(x = Age, y = y, label = lab),
              hjust = 0, nudge_x = 0.3, size = label_size, fontface = "bold",
              inherit.aes = FALSE, colour = pal_label[as.character(band_labels$Band)]) +
    facet_wrap(~Scenario, labeller = strip_fn) +
    scale_fill_manual(values = pal_outcome) +
    x_shared +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0.035, 0.02))) +
    labs(title = bottom_title, subtitle = bottom_subtitle, caption = bottom_caption,
         y = "Share of the\nstarting cohort") +
    base +
    theme(strip.text   = element_text(size = base_size * 0.62, face = "plain",
                                      colour = "grey35", hjust = 0),
          plot.caption = element_text(hjust = 0, colour = "grey35",
                                      size = base_size * 0.62, margin = margin(t = 10)))

  (p_top / p_bottom) + plot_layout(heights = heights)
}
