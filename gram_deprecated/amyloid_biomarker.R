# AMYLOID BETA BIOMARKER AND ALZHEIMER'S DISEASE


###### Set up ######
# It seems from glancing through Men's data extraction sheet that studies more often have the PPA/NPA value rather than sens/spec.
# I use the following optimization chunk to calculate specificity and sensitivity values for known prevalence and target PPA and NPA values.
# Define the objective function
objective_function <- function(params) {
  sensitivity <- params[1]
  specificity <- params[2]
  
  # Calculate PPA and NPA
  ppa <- (sensitivity * prevalence) / ((sensitivity * prevalence) + ((1 - specificity) * (1 - prevalence)))
  npa <- (specificity * (1 - prevalence)) / ((specificity * (1 - prevalence)) + ((1 - sensitivity) * prevalence))
  
  # Calculate the absolute differences from desired values
  ppa_diff <- abs(ppa - 0.77)
  npa_diff <- abs(npa - 0.78)
  
  # Return the sum of absolute differences (to minimize)
  return(ppa_diff + npa_diff)
}

# Initial guess for sensitivity and specificity
initial_guess <- c(0.5, 0.5)

# Run optimization
result <- optim(initial_guess, objective_function)

# Extract optimized sensitivity and specificity
optimized_params <- result$par
optimized_sensitivity <- optimized_params[1]
optimized_specificity <- optimized_params[2]

# Calculate optimized PPA and NPA
optimized_ppa <- (optimized_sensitivity * prevalence) / ((optimized_sensitivity * prevalence) + ((1 - optimized_specificity) * (1-prevalence)))
optimized_npa <- (optimized_specificity * (1 - prevalence)) / ((optimized_specificity * (1 - prevalence)) + ((1 - optimized_sensitivity) * prevalence))

# Print results
cat("Optimized Sensitivity:", optimized_sensitivity, "\n")
cat("Optimized Specificity:", optimized_specificity, "\n")
cat("Optimized PPA:", optimized_ppa, "\n")
cat("Optimized NPA:", optimized_npa, "\n")




###### Test model #####


# Start by defining a test population 
## sex 0 = male, 1 = female
## ab (amyloid-beta) 0 = less than cutoff, 1 = equal to or more than cutoff
df <- data.frame(id = c(1,2), sex = c(0,1), ab_test = c(0,1), ab_true = NA, alz_true = NA)

# True prevalence of amyloid plaques (NOT ALZHEIMERS) in population - based on gold standard (in this case, PET)
prevalence <- 0.375

# Test performance for plasma amyloid-beta test, compared to gold standard (at some point can incorporate different cut-offs)
# This is derived in the "Set up" chunk above - this parameterization could be built into the active model, or we could change for different tests etc if already available.
sensitivity <- 0.58
specificity <- 0.90
ppa <- (sensitivity * prevalence) / ((sensitivity * prevalence) + ((1-specificity) * (1-prevalence)))
npa <- (specificity * (1 - prevalence)) / ((specificity * (1 - prevalence)) + ((1 - sensitivity) * prevalence))

# Determine true amyloid status of agents (i.e., people) based on plasma biomarker test result
df$ab_true[df$ab_test == 1] <- rbinom(sum(df$ab_test), size = 1, prob = ppa)
df$ab_true[df$ab_test == 0] <- 1 - rbinom(sum(1-df$ab_test), size = 1, prob = npa)

# If they have the plasma biomarker, what is the probability that they currently have or will develop Alzheimer's Disease?
# For now, assume 10% -- this will need to be a formula accounting for other conditions that might cause amyloid in plasma
# QUESTION: Is it true that age is a factor probably of being amyloid positive, but NOT for the subsequent conditional probability of amyloid-->Alzheimers?
prob_alz_if_ab <- 0.1
df$alz_true <- df$ab_true * rbinom(n = nrow(df), size = 1, prob = prob_alz_if_ab)
df


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8573146/: Study has regression coefficients for amyloid beta concentration and results of neuropsych tests. 
# If we want to use something like this, we will need individual-level data on CONTINUOUS amyloid beta plasma concentration.
# That would be ideal to also model the increase in concentration over time, rather than relying on a binary variable.
# How likely is it that we could get individual-level data from a longitdunak AD cohort?