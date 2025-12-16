make_gram_pie_chart <- function(output_array, cycles = 16:31, plot_title = "Cumulative Testing Status, Ages 65-79") {

  # Subset output array to only include
  # a. alive at 65
  # b. never identified outside of BHA pathway (DX never 1)
  # c. attributes necessary for calculations
  
  alive_at_first_test <- output_array[cycles[1], "ALIVE" ,] == 1
  has_hcare <- colSums(output_array[cycles, "HCARE", ], na.rm = TRUE) > 0
  never_dx <- colSums(output_array[cycles-1, "DX", ], na.rm = TRUE) == 0
  select_ids <- alive_at_first_test & has_hcare & never_dx  # 64,024 people left in u1bhapos
  
  output_subset <- output_array[ , , select_ids]
  
  bha_mat <- output_array[cycles, "BHA", select_ids]
  pcp_mat <- output_array[cycles, "PCP", select_ids]
  
  syn_mat <- output_array[cycles, "SYN", select_ids]
  
  
  # Combine BHA and PCP into single outcome matrix
  outcome_mat <- pcp_mat
  outcome_mat[is.na(pcp_mat)] <- bha_mat[is.na(pcp_mat)]
  
  # Calculate outcome at each test (contemporaneous)
  cont_outcome_mat <- matrix(NA, nrow = cycles, ncol = sum(select_ids))
  
  cont_outcome_mat[outcome_mat == 1 & syn_mat == 1] <- "TP"
  cont_outcome_mat[outcome_mat == 0 & syn_mat == 0] <- "TN"
  cont_outcome_mat[outcome_mat == 0 & syn_mat == 1] <- "FN"
  cont_outcome_mat[outcome_mat == 1 & syn_mat == 0] <- "FP"
  
  cont_outcome_mat[outcome_mat < 0 & syn_mat == 0] <- "TN*"  # not tested and healthy -- default true negative
  cont_outcome_mat[outcome_mat < 0 & syn_mat == 1] <- "FN*"  # not tested and impaired -- default false negative
  
  
  # Calculate aggregate outcomes
  ever_impaired <- colSums(syn_mat == 1, na.rm = TRUE) > 0
  ever_positive <- colSums(outcome_mat == 1, na.rm = TRUE) > 0
  never_tested <- colSums(bha_mat == 1 | bha_mat == 0, na.rm = TRUE) == 0
  
  agg_outcome_vec <- rep(NA, times = sum(select_ids))
  
  agg_outcome_vec[ever_impaired  & ever_positive] <- "TP"
  agg_outcome_vec[!ever_impaired & !ever_positive] <- "TN"
  agg_outcome_vec[ever_impaired  & !ever_positive] <- "FN"
  agg_outcome_vec[!ever_impaired & ever_positive] <- "FP"
  
  agg_outcome_vec[!ever_impaired & never_tested] <- "TN*"  # healthy and not tested -- default true negative
  agg_outcome_vec[ever_impaired  & never_tested] <- "FN*"   # impaired and not tested -- default false negative
  
  counts <- c(
    "True Positive" = sum(agg_outcome_vec == "TP"),
    "False Positive" = sum(agg_outcome_vec == "FP"),
    "True Negative" = sum(agg_outcome_vec == "TN"),
    "False Negative" = sum(agg_outcome_vec == "FN"),
    "Never tested & Healthy" = sum(agg_outcome_vec == "TN*"),
    "Never tested & Impaired" = sum(agg_outcome_vec == "FN*")
  )
  
  bha_df_pies <- data.frame(status = names(counts), count = as.numeric(counts))
  
  piechart <- ggplot(bha_df_pies, aes(x = "", y = count, fill = status)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    labs(title = plot_title, fill = "") +
    theme_void() +
    scale_fill_manual(values = c("red", "pink", "skyblue", "orange", "forestgreen", "goldenrod"))
  
  # test_trajectory <- data.frame(
  #   consistently_tp_tn = consistently_tp_tn,
  #   fp_to_tn = fp_to_tn,
  #   fn_to_tp = fn_to_tp
  # )
  return(list(df = bha_df_pies, piechart = piechart)) 
              # test_trajectory_stats = test_trajectory))
}

