# Tau-PET based Braak progression model
# v1.0.0

# Load libraries
library(ggplot2)
library(reshape2)

# Read in transition probability data (see braak_progression_rates.R and tau_histology_calibration.R for estimation of these values)
trans_probs_abbr <- readRDS("tau_tangle_trans_probs_calibrated.RDS")

braak_stages <- c(0, 1, 2, 3, 4, 5, 6)       # 6 Braak stages, sequential, plus one Braak 0 stage
braak_stage_distribution <- c(0.09, 0.08, 0.10, 0.11, 0.12, 0.20, 0.30)  # Distribution assumed based on most recent NACC neuropathology data
pet_sensitivities <- c(1, 0.5, 0.5, 0.8, 0.8, 0.99, 1) # Sensitivity of PET based on true Braak stage - assumed (arbitrary placeholders)


# Custom function that depicts change in PET-based Braak stage over time, with n_i individuals over n_t years
tau_model <- function(n_i, n_t, seed=20240522) {

  
  # Set a random seed
  set.seed(seed)
  
  # Create a data frame with n_i individuals to hold individual-level attributes (in this version )
  agents <- data.frame(
    id = seq(1:n_i),                 # Assign each agent and ID
    sex = rbinom(n_i, 1, 0.5),       # 0 = male, 1 = female, 50/50
    educ = rbinom(n_i, 1, 0.6),      # 0 = less than college (40%), 1 = college graduate (60%)
    age =  as.integer(rnorm(n_i, mean=60, sd=10)),
    ad_path = rbinom(n_i, 1, 0.75)   # 0 = no AD pathology, 1 = has AD pathology, arbitrarily set to 75%   
  )
  
  # Assign initial Braak stages
  n_ad_path <- sum(agents$ad_path)                                 # Count number of people with AD pathology
  agents$true_braak[agents$ad_path == 0] <- braak_stages[1]        # If not on AD path, assign B0
  agents$true_braak[agents$ad_path == 1] <- sample(braak_stages, size=n_ad_path, replace = TRUE, prob = braak_stage_distribution)  # If on AD path, randomly assign a Braak stage
  
  # Calculate initial PET-based Braak stages
  accurate_pet_flag <- runif(n_i) <= pet_sensitivities[match(agents$true_braak, braak_stages)]
  agents$pet_braak <- braak_stages[match(agents$true_braak, braak_stages) - (1-accurate_pet_flag)]
  
  
  # Create a 3D array to hold model data. This is a collection of matrices, with each matrix representing a "slice" in time.
  # Each slice has individuals in rows and different outcomes, with different time-points stored in the third dimension. 
  # NB for now we don't use any of the downstream outcomes
  variables <- c("id", "sex", "educ", 
                 "age", "ad_path", 
                 "true_braak", "pet_braak", 
                 "mmse", "ab42", "atau", 
                 "dementia", "qol", "cost")
  model <- array(NA, c(n_i, length(variables),  n_t+1),
                 dimnames = list(paste0("ind ", 1:n_i), variables, paste0("cyc ", 0:n_t)))
  
  # Assign starting (cycle 0) values from the agents data frame created above
  # Note that ID, Sex, Educ, and AD path are time-independent, therefore they're written onto all t matrices
  model[ ,"id" , ] <- agents$id
  model[ ,"sex" , ] <- agents$sex
  model[ ,"educ" , ] <- agents$educ
  model[ ,"age" , 1] <- agents$age
  model[ ,"ad_path" , ] <- agents$ad_path
  model[ ,"true_braak" , 1] <- agents$true_braak
  model[ ,"pet_braak" , 1] <- agents$pet_braak
  
  # Run through time cycles (annual)
  for (t in 1:n_t) {
    
  # Calculate next-cycle Braak stage. 
  # Note, Braak stages are thought to monotonically and sequentially increase. 
  # So only two probabilities exist at each cycle: either maintain previous stage, or progress to the next stage
    
  # First, generate vector of length n_i with probabilities of staying in the stage Braak stage
  ##### If not on AD path, prob=1 for staying B0. If on AD path, the appropriate probability for the relevant Braak stage is looked up from trans_probs_abbr
  # Next, decide if they will stay or progress based on these probabilities with a random dice roll
  ##### If the random number is smaller than the probability of moving, they stay. Otherwise, they move one stage up. Store these outcomes in progress_flag.
  # Lastly, get the index for current Braak stage and add progress_flag (0 if staying, 1 if progressing). Get new Braak stage for t+1.
  prob_no_change <- ifelse(model[ ,"ad_path", t] == 0, 1,
                           trans_probs_abbr$probability[match(model[ , "true_braak", t], trans_probs_abbr$current_braak)])
  progress_flag <- runif(n_i) > prob_no_change
  model[ , "true_braak", t+1] <- braak_stages[match(model[ , "true_braak", t], braak_stages) + progress_flag]
  
  # Find PET-based Braak stage (as done above for agents$pet_braak, but using updated true Braak staging)
  accurate_pet_flag <- runif(n_i) <= pet_sensitivities[match(model[ ,"true_braak", t+1], braak_stages)]
  model[ , "pet_braak", t+1] <- braak_stages[match(model[ , "true_braak", t+1], braak_stages) - (1-accurate_pet_flag)]

  } 
  
  # Create simulation trace by reshaping model data
  # aggregated_data is the trace of proportion of individuals at each stage in each cycle in table format)
  # trace is long version of the same, necessary for ggplot below
  aggregated_true_braak <- prop.table(table(melt(model[ , "true_braak", ])$Var2, melt(model[ , "true_braak", ])$value), margin = 1)
  trace_true_braak <- melt(aggregated_true_braak, value.name = "Proportion")  
    
  # Plot stacked area chart of proportion of each Braak stage over time
  plot_true_tau <- ggplot(trace_true_braak, aes(x = Var1, y = Proportion, group = Var2, fill = Var2)) +
    geom_area() +
    labs(title = "Distribution of True Braak Stages Over Time", x = "Time", y = "Proportion of Individuals", fill = "True Braak Stage") +
    theme_minimal()
  
  
  aggregated_pet_braak <- prop.table(table(melt(model[ , "pet_braak", ])$Var2, melt(model[ , "pet_braak", ])$value), margin = 1)
  trace_pet_braak <- melt(aggregated_pet_braak, value.name = "Proportion")  
  
  # Plot stacked area chart of proportion of each Braak stage over time
  plot_pet_tau <- ggplot(trace_pet_braak, aes(x = Var1, y = Proportion, group = Var2, fill = Var2)) +
    geom_area() +
    labs(title = "Distribution of PET-based Braak Stages Over Time", x = "Time", y = "Proportion of Individuals", fill = "PET-based Braak Stage") +
    theme_minimal()
  
  # Return model results including the full model data (3D array), the stacked area plot, and the simulation trace table
  return(list(model_output = model, model_plots = list(plot_true_tau, plot_pet_tau), model_traces = list(aggregated_true_braak, aggregated_pet_braak)))
  
}


model1 <- tau_model(n_i=1000, n_t=10)

model1$model_plots[1]
model1$model_plots[2]
model1$model_output
true_braak <- dcast(as.data.frame(model1$model_traces[1]), Var1 ~ Var2, value.var = "Freq") 
pet_braak <- dcast(as.data.frame(model1$model_traces[2]), Var1 ~ Var2, value.var = "Freq")


# Calculate modeled two-year probabilities
sum_probs_stay <- vector(mode = "numeric", length = 7)

for (row in (1:(nrow(pet_braak)-2))) {
  sum_probs_stay <- sum_probs_stay + pmin(1, (pet_braak[(row+2), -1] / pet_braak[row, -1]))
}
sum_probs_stay / (nrow(pet_braak)-2)

# Calculate annual probabilities
sum_probs_stay <- vector(mode = "numeric", length = 7)

for (row in (1:(nrow(pet_braak)-2))) {
  sum_probs_stay <- sum_probs_stay + pmin(1, (pet_braak[(row+1), -1] / pet_braak[row, -1]))
}
sum_probs_stay / (nrow(pet_braak)-1)


# Next step ---
# run the calibration on these simulated two-year probs to get more accurate results!
# this requires reformatting of the outputs similar to what I did in the calibration script
# Also for the future -- the model output is an array, which only holds one type of element (in this case, character).
# This needs to be rearranged in the future so that we can utilize multiple types of arguments (alternatively could constantly convert on the fly as needed?)