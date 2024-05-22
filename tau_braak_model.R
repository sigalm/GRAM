# Tau-PET based Braak progression model
# v1.0.0

# Load libraries
library(ggplot2)
library(reshape2)

# Read in transition probability data (see braak_progression_rates.R for estimation of these values)
trans_probs_abbr <- readRDS("trans_probs_abbr.rds")

# Custom function that depicts change in PET-based Braak stage over time, with n_i individuals over n_t years
tau_model <- function(n_i, n_t, seed=20240522) {

  
  # Set a random seed
  set.seed(seed)
  
  # Create a data frame with n_i individuals to hold individual-level attributes (in this version )
  agents <- data.frame(
    id = seq(1:n_i),                 # Assign each agent and ID
    sex = rbinom(n_i, 1, 0.5),       # 0 = male, 1 = female, 50/50
    educ = rbinom(n_i, 1, 0.6),      # 0 = less than college (40%), 1 = college graduate (60%)
    age =  as.integer(rnorm(n_i, 60, sd=10)),
    ad_path = rbinom(n_i, 1, 0.75)   # 0 = no AD pathology, 1 = has AD pathology, arbitrarily set to 50%   
  )
  
  

  # Assign initial Braak stages
  n_ad_path <- sum(agents$ad_path)                            # Count number of people with AD pathology
  braak_stages <- c("B0","B1","B2","B3","B4","B5","B6")       # 6 Braak stages, sequential, plus one no-Braak stage
  agents$braak[agents$ad_path == 0] <- braak_stages[1]        # If not on AD path, assign B0
  agents$braak[agents$ad_path == 1] <- sample(braak_stages, size=n_ad_path, replace = TRUE, prob = rep(x=1/7, 7))  # If on AD path, randomly assign a Braak stage (for now, assume uniform equal distribution) 
  
  
  # Create a 3D array to hold model data. This is a collection of matrices, with each matrix representing a "slice" in time.
  # Each slice has individuals in rows and different outcomes, with different time-points stored in the third dimension. 
  # NB for now we don't use any of the downstream outcomes
  variables <- c("id", "sex", "educ", 
                 "age", "ad_path", 
                 "braak", "mmse", "ab42", "atau", 
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
  model[ ,"braak" , 1] <- agents$braak
  
  
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
                           trans_probs_abbr$probability[match(model[ , "braak", t], trans_probs_abbr$baseline)])
  progress_flag <- runif(n_i) > prob_no_change
  model[ , "braak", t+1] <- braak_stages[match(model[ , "braak", t], braak_stages) + progress_flag]

  } 
  
  # Create simulation trace by reshaping model data
  # aggregated_data is the trace of proportion of individuals at each stage in each cycle in table format)
  # trace is long version of the same, necessary for ggplot below
  aggregated_data <- prop.table(table(melt(model[ , "braak", ])$Var2, melt(model[ , "braak", ])$value), margin = 1)
  trace <- melt(aggregated_data, value.name = "Proportion")  
    
  # Plot stacked area chart of proportion of each Braak stage over time
  plot <- ggplot(trace, aes(x = Var1, y = Proportion, group = Var2, fill = Var2)) +
    geom_area() +
    labs(title = "Distribution of Braak Stages Over Time", x = "Time", y = "Proportion of Individuals", fill = "Braak Stage") +
    theme_minimal()
  print(plot)
  
  # Return model results including the full model data (3D array), the stacked area plot, and the simulation trace table
  return(list(model_output = model, model_plot = plot, model_trace = aggregated_data))
  
}


model1 <- tau_model(n_i=1000,n_t=20)

model1$model_plot
model1$model_output
model1$model_trace
