
# Braak annual progression rates optimization
# Data from: Therriault 2022, Fig 1c  (https://www.nature.com/articles/s43587-022-00204-0#data-availability)

# I define the possible Braak stages, then assign the one-year transition probabilities from the table on the left-hand side of figure.
# Using these 1-year probabilities, I calculate the two-year probabilities in a function
# Then, I run the function and the initial one-year parameters through an optimization function.
## This allows me to calibrate the one-year prob values that get me closest to the two-year progression rates. 
# I then save the "optimized" (i.e., calibrated) one-year parameters into a table.
# Note that at each stage, only 2 transitions are possible: 
## either proceed one step, or stay the same (authors note that any "improvement" in stage is likely a misclassification, so I included those with the "stay put" group)
# These transitions are the complements of one another, so we can simplify the transition probabilities table to only include the "stay put" probabilities at each stage
# The prob of progressing one step will be the complement of these probabilities. 

braak_stages <- c("B0","B1","B2","B3","B4","B5","B6")

b00 <- 0.956
b11 <- 1
b22 <- 0.767
b33 <- 0.25
b44 <- 0.833
b55 <- 1
b66 <- 1

# Define a function to calculate two-year transition probabilities
calc_two_year_probs <- function(params) {
  b00 <- params[1]
  b01 <- 1 - b00
  b11 <- params[2]
  b12 <- 1 - b11
  b22 <- params[3]
  b23 <- 1 - b22
  b33 <- params[4]
  b34 <- 1 - b33
  b44 <- params[5]
  b45 <- 1 - b44
  b55 <- params[6]
  b56 <- 1 - b55
  b66 <- params[7]
  
  # Create the annual transition matrix
  annual_matrix <- matrix(0, nrow = 7, ncol = 7)
  annual_matrix[1,] <- c(b00, b01, 0, 0, 0, 0, 0)
  annual_matrix[2,] <- c(0, b11, b12, 0, 0, 0, 0)
  annual_matrix[3,] <- c(0, 0, b22, b23, 0, 0, 0)
  annual_matrix[4,] <- c(0, 0, 0, b33, b34, 0, 0)
  annual_matrix[5,] <- c(0, 0, 0, 0, b44, b45, 0)
  annual_matrix[6,] <- c(0, 0, 0, 0, 0, b55, b56)
  annual_matrix[7,] <- c(0, 0, 0, 0, 0, 0, b66)
  
  # Calculate the two-year transition matrix
  two_year_matrix <- annual_matrix %*% annual_matrix
  return(two_year_matrix)
}

# Define the given two-year transition probabilities
# (This should be replaced with your actual two-year transition probabilities)
given_two_year_probs <- matrix(c(
  0.934, 0.033, 0.022, 0, 0, 0, 0,
  0, 0.66, 0.33, 0, 0, 0, 0,
  0, 0, 0.690, 0.230, 0.080, 0, 0,
  0, 0, 0, 0, 1, 0, 0,
  0, 0, 0, 0, 0.800, 0.200, 0,
  0, 0, 0, 0, 0, 0.500, 0.500,
  0, 0, 0, 0, 0, 0, 1
), nrow = 7, byrow = TRUE)

# Define the objective function
objective_function <- function(params) {
  two_year_matrix <- calc_two_year_probs(params)
  sum((two_year_matrix - given_two_year_probs)^2)
}

initial_params <- c(b00, b11, b22, b33, b44, b55, b66)

opt_result <- optim(par = initial_params, fn = objective_function, method = "L-BFGS-B",
                    lower = rep(0, 7), upper = rep(1, 7))

# Extract optimized parameters
optimized_params <- opt_result$par
optimized_params

optimized_matrix <- matrix(0, nrow = 7, ncol = 7)
optimized_matrix[1,] <- c(optimized_params[1], 1 - optimized_params[1], 0, 0, 0, 0, 0)
optimized_matrix[2,] <- c(0, optimized_params[2], 1 - optimized_params[2], 0, 0, 0, 0)
optimized_matrix[3,] <- c(0, 0, optimized_params[3], 1 - optimized_params[3], 0, 0, 0)
optimized_matrix[4,] <- c(0, 0, 0, optimized_params[4], 1 - optimized_params[4], 0, 0)
optimized_matrix[5,] <- c(0, 0, 0, 0, optimized_params[5], 1 - optimized_params[5], 0)
optimized_matrix[6,] <- c(0, 0, 0, 0, 0, optimized_params[6], 1 - optimized_params[6])
optimized_matrix[7,] <- c(0, 0, 0, 0, 0, 0, optimized_params[7])

rownames(optimized_matrix) <- colnames(optimized_matrix) <- braak_stages
print(optimized_matrix)
rowSums(optimized_matrix)

# Melt into dataframe for easier application
library(reshape2)
trans_probs <- melt(optimized_matrix, varnames = c("baseline","followup"),value.name = "probability")
trans_probs <- trans_probs[order(trans_probs$baseline, trans_probs$followup), ]
rownames(trans_probs) <- NULL

trans_probs_abbr <- trans_probs[trans_probs$baseline == trans_probs$followup, ]



saveRDS(trans_probs_abbr, "trans_probs_abbr.rds")
