# Annual changes in underlying Braak histopathology

# Assumption 1: Hierarchical progression (i.e., cannot skip stages)
# Assumption 2: Monotonically worsening progression -- no improvement in Braak stage

# Create vector of Braak stages
braak_stages <- c(0,1,2,3,4,5,6)        # 6 Braak stages, sequential, plus one Braak 0 stage

# Define starting distribution in the population
braak_stage_distribution <- c(0.09, 0.08, 0.10, 0.11, 0.12, 0.20, 0.30)  # Distribution assumed based on most recent NACC neuropathology data

# Define sensitivity of PET for each Braak stage. For now, we disregard specificity (i.e., assume no false-positive staging)
# Vector elements follow the same order as braak_stages
pet_sensitivities <- c(1, 0.5, 0.5, 0.8, 0.8, 0.99, 1)

# Write a function to calculate annual and two-year probabilities of staying in the same Braak stage
#    given unobserved underlying histological probabilities and known initial distribution and test performance
pet_transitions <- function(probs_stay) {
  
  probs_stay <- c(probs_stay, 1)
  probs_progress <- 1 - probs_stay
  
  # Calculate true histopathology over time
  true_hist_progression <- data.frame(T0 = braak_stage_distribution,
                                      T1 = NA,
                                      T2 = NA,
                                      row.names = braak_stages)
  
  for (i in seq(braak_stages)) {
    if(i==1) {
      true_hist_progression$T1[i] <- true_hist_progression$T0[i] * probs_stay[i]
      true_hist_progression$T2[i] <- true_hist_progression$T1[i] * probs_stay[i]
    } else {
      true_hist_progression$T1[i] <- true_hist_progression$T0[i] * probs_stay[i] + true_hist_progression$T0[i-1] * probs_progress[i-1]
      true_hist_progression$T2[i] <- true_hist_progression$T1[i] * probs_stay[i] + true_hist_progression$T1[i-1] * probs_progress[i-1]
    }
  }
  
  pet_accurate <- data.frame(T0 = true_hist_progression$T0 * pet_sensitivities,
                             T1 = true_hist_progression$T1 * pet_sensitivities,
                             T2 = true_hist_progression$T2 * pet_sensitivities)
  
  # Since PET less likely to detect earlier stages, assume that inaccurately staged PET results are misclassified by 1 level
  # i.e., if 0.08 is B1 with 50% sensitivity, 0.04 will be accurately staged as B1, the rest will be B0
  # if 0.10 is B2 with 60% sensitivity, 0.06 will be accurately staged as B2, the rest will be B1 (NOT B0)
  pet_downstaged <- data.frame(T0 = true_hist_progression$T0 * (1-pet_sensitivities),
                               T1 = true_hist_progression$T1 * (1-pet_sensitivities),
                               T2 = true_hist_progression$T2 * (1-pet_sensitivities))
  
  pet_overall <- pet_accurate + rbind(pet_downstaged[2:7, ], c(0,0,0))
  
  pet_annual_probs_1 <- pmin(1, pet_overall$T1 / pet_overall$T0)
  pet_annual_probs_2 <- pmin(1, pet_overall$T2 / pet_overall$T1)
  
  pet_annual_probs_mean <- (pet_annual_probs_1 + pet_annual_probs_2) / 2
  
  pet_two_year_probs <- pmin(1, pet_overall$T2 / pet_overall$T0)
  
  result <- data.frame(annual_probs = pet_annual_probs_mean,
                       two_year_probs = pet_two_year_probs,
                       row.names = braak_stages)
  
  return(result)
}

# Define target "stay put" probabilities
target_data <- data.frame(annual = c(0.973, 0.807, 0.833, 0.151, 0.913, 0.709),
                          two_year = c(0.934, 0.660, 0.690, 0, 0.80, 0.5))
                 
# Define a function to calculate discrepancy between predictions and target data (in this case, sum of squared errors [SSE])
objective_function <- function(probs_stay) {
  predicted_probs <- pet_transitions(probs_stay)
  sse <- sum((predicted_probs$two_year_probs[-7] - target_data$two_year)^2)
  #sse_two_year <- sum((predicted_probs$two_year_probs - target_data$two_year_probs)^2)
  #total_sse <- sse_annual + sse_two_year
}

# Define initial guess for the annual probabilities of staying in the same stage
# First guess: All stages have 0.8 probability of staying the same next year, except for B6 which is 1 by definition
initial_probs <- c(0.9, 0.5, 0.5, 0.5, 0.5, 0.5) 

# Perform optimization
# Lower and upper arguments are vectors holding the lower and upper pounds of the parameters we're varying (in this case all bound at 0-1 since probabilities)
# The method is...
opt_results <- optim(par = initial_probs,
                     fn = objective_function,
                     method = "L-BFGS-B",
                     lower = rep(0, 6), upper = rep(1, 6))

calibrated_probs <- opt_results$par
print(calibrated_probs)
calibrated_results <- pet_transitions(calibrated_probs)

# Investigate residuals to evaluate calibration success
residual_error <- opt_results$value
residuals <- calibrated_results$two_year_probs[-7] - target_data$two_year
residuals_df <- data.frame(Braak_Stage = braak_stages[-7], Residuals = residuals)

library(ggplot2)
ggplot(residuals_df, aes(x = Braak_Stage, y = Residuals)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Residuals of Two-Year Probabilities by Braak Stage",
       x = "Braak Stage",
       y = "Residuals")


# If accepting calibration, save transition probabilities in an RDS file
probs_out <- data.frame(current_braak = braak_stages, probability = calibrated_results$annual_probs)
saveRDS(probs_out, file = "tau_tangle_trans_probs_calibrated.RDS")
