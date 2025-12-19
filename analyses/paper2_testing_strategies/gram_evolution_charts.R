library(ggplot2)
library(tidyr)











add_early_dx <- function(output_array, dx_var) {
  clinical_dx <- output_array[,"DX",]
  scenario_tp_direct <- output_array[,"TP_direct",]
  scenario_tp_indirect <- output_array[,"TP_indirect",]
  scenario_fp <- output_array[,"FP",]
  syn_status <- output_array[,"SYN",]
  
  # Helper function to find the first occurrence time for each person (column)
  get_first_occurrence_time <- function(event_matrix) {
    apply(event_matrix, 2, function(person_vector) {
      time <- min(which(person_vector == 1))
      ifelse(is.infinite(time), NA, time)
    })
  }
  
  # --- Pre-calculate key event times and flags for each person ---
  
  # Time of first clinical diagnosis (DX)
  time_clinical_dx <- get_first_occurrence_time(clinical_dx)
  ever_clinical_dx <- !is.na(time_clinical_dx)
  
  # Time of first scenario True Positive (direct or indirect)
  scenario_tp <- scenario_tp_direct | scenario_tp_indirect
  time_scenario_tp <- get_first_occurrence_time(scenario_tp)
  
  # Time of first scenario positive test (TP direct or FP)
  scenario_pos_test <- scenario_tp_direct | scenario_fp
  time_pos_test <- get_first_occurrence_time(scenario_pos_test)
  
  # Time of first False Positive (FP)
  time_fp <- get_first_occurrence_time(scenario_fp)
  
  # --- Identify individuals based on the 5 conditions ---
  
  # Condition 1 & 2: Gained TP in scenario, but would never get clinical DX
  # This is a per-person attribute, not a time-varying one.
  # We check if they have a scenario TP and would NOT have a clinical DX.
  gained_tp_no_clinical_dx <- !ever_clinical_dx & !is.na(time_scenario_tp)
  
  # Identify which type of TP it was (direct or indirect) at the time of the event.
  time_tp_direct <- get_first_occurrence_time(scenario_tp_direct)
  time_tp_indirect <- get_first_occurrence_time(scenario_tp_indirect)
  
  # early_dx_1: Event happens at the time of TP_direct for those who never get a clinical DX.
  early_dx_1 <- matrix(FALSE, nrow = nrow(clinical_dx), ncol = ncol(clinical_dx))
  for (person_idx in which(gained_tp_no_clinical_dx & !is.na(time_tp_direct))) {
    early_dx_1[time_tp_direct[person_idx], person_idx] <- TRUE
  }
  
  # early_dx_2: Event happens at the time of FP for those who get a TP_indirect but never a clinical DX.
  early_dx_2 <- matrix(FALSE, nrow = nrow(clinical_dx), ncol = ncol(clinical_dx))
  # The condition identifies people who had a TP_indirect. The event is recorded at the time of the FP.
  for (person_idx in which(gained_tp_no_clinical_dx & !is.na(time_tp_indirect))) {
    # Ensure the person had an FP, then record the event at that time.
    if (!is.na(time_fp[person_idx])) {
      early_dx_2[time_fp[person_idx], person_idx] <- TRUE
    }
  }
  
  # Condition 3 & 4: Gained TP in scenario earlier than clinical DX
  # This applies only to people who get both a scenario TP and a clinical DX.
  gained_tp_earlier <- ever_clinical_dx & !is.na(time_scenario_tp) & (time_scenario_tp < time_clinical_dx)
  
  # early_dx_3: Event at time of TP_direct if it's earlier than clinical DX.
  early_dx_3 <- matrix(FALSE, nrow = nrow(clinical_dx), ncol = ncol(clinical_dx))
  for (person_idx in which(gained_tp_earlier & !is.na(time_tp_direct) & (time_tp_direct < time_clinical_dx[person_idx]))) {
    early_dx_3[time_tp_direct[person_idx], person_idx] <- TRUE
  }
  
  # early_dx_4: Event at time of FP if it leads to a TP_indirect that is earlier than clinical DX.
  early_dx_4 <- matrix(FALSE, nrow = nrow(clinical_dx), ncol = ncol(clinical_dx))
  # The condition identifies people whose TP_indirect was earlier than their clinical DX. The event is recorded at the time of the FP.
  for (person_idx in which(gained_tp_earlier & !is.na(time_tp_indirect) & (time_tp_indirect < time_clinical_dx[person_idx]))) {
    # Ensure the person had an FP, then record the event at that time.
    if (!is.na(time_fp[person_idx])) {
      early_dx_4[time_fp[person_idx], person_idx] <- TRUE
    }
  }
  
  # Condition 5: Received clinical DX before any scenario positive test (TP direct or FP)
  # This can happen if the person has a clinical DX but no positive test, or if the DX comes first.
  clinical_dx_before_pos_test <- ever_clinical_dx & (is.na(time_pos_test) | (time_clinical_dx < time_pos_test))
  
  # early_dx_5: Event happens at the time of the clinical DX.
  early_dx_5 <- matrix(FALSE, nrow = nrow(clinical_dx), ncol = ncol(clinical_dx))
  for (person_idx in which(clinical_dx_before_pos_test)) {
    early_dx_5[time_clinical_dx[person_idx], person_idx] <- TRUE
  }
  
  # --- Append new variables to the output array ---
  output_modified <- abind(output_array, 
                           "early_dx_1" = early_dx_1, 
                           "early_dx_2" = early_dx_2,
                           "early_dx_3" = early_dx_3,
                           "early_dx_4" = early_dx_4,
                           "early_dx_5" = early_dx_5,
                           along = 2)
  
  return(output_modified)
}



