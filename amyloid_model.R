
my_model <- function(n_i, n_t, seed=12345) {

set.seed(seed)

agents <- data.frame(
  id = seq(1:n_i),                 # Assign each agent and ID
  sex = rbinom(n_i, 1, 0.5),       # 0 = male, 1 = female, 50/50
  educ = rbinom(n_i, 1, 0.6),      # 0 = less than college (40%), 1 = college graduate (60%)
  age =  as.integer(rnorm(n_i, 60, sd=10)),
  ad_path = rbinom(n_i, 1, 0.35)   # 0 = no AD pathology (65%), 1 = has AD pathology (35%))   
)

## Assign dependent attributes

n_ad_path <- sum(agents$ad_path)
agents$plaq_density <- NA
agents$plaq_density[agents$ad_path == 1] <- rnorm(n_ad_path, mean = 0.78, sd = 0.02)
agents$plaq_density[agents$ad_path == 0] <- rnorm(nrow(agents) - n_ad_path, mean = 0.74, sd = 0.03)


# plaque density and rate of change data from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8055327/
# They stratify the baseline amyloid density and annual rate of change by subsequent conversion status to amyloid beta positive.
# Since I'm not clear on what constitutes AD pathology, I assumed converters to be AD_path positive, and non-converters to be AD_path negative.


##### Other vars to be defined later #####
# dementia <- to be defined as MMSE below a certain point? 
# qol <- 
# cost <- 


## Assign non-updated attributes

# mmse <- from Geena model
# amyloid_b42 <- need to find some data

#####

# Create data frames to hold outcomes
variables <- c("id", "sex", "educ", 
               "age", "ad_path", 
               "plaq_density", "mmse", "ab42", "atau", 
               "dementia", "qol", "cost")
model <- array(NA, c(n_i, length(variables),  n_t+1),
               dimnames = list(paste0("ind ", 1:n_i), variables, paste0("cyc ", 0:n_t)))

model[ ,"id" , ] <- agents$id
model[ ,"sex" , ] <- agents$sex
model[ ,"educ" , ] <- agents$educ
model[ ,"age" , 1] <- agents$age
model[ ,"ad_path" , ] <- agents$ad_path
model[ ,"plaq_density" , 1] <- agents$plaq_density


# Annual rate of amyloid beta deposition (via PET):
# In non-converters: 0 +- 0.005
# In converters: 0.014 +- 0.004

mean_amyloid_annual_change <- 0.014
sd_amyloid_annual_change <- 0.004

mean_amyloid_annual_change_0 <- 0
sd_amyloid_annual_change0 <- 0.005

for (t in 1:n_t) {
  # If has AD pathology, use the annual change for converters
  model[model[ , "ad_path", t] == 1 , "plaq_density", t+1] <-  
    model[model[ , "ad_path", t] == 1, "plaq_density", t] * (1 + rnorm(n_ad_path, mean_amyloid_annual_change, sd_amyloid_annual_change))
  
  # If does not have AD pathology, use the annual change for non-converters
  model[model[ , "ad_path", t] == 0 , "plaq_density", t+1] <-  
    model[model[ , "ad_path", t] == 0, "plaq_density", t] * (1 + rnorm(n_i-n_ad_path, mean_amyloid_annual_change_0, sd_amyloid_annual_change0))
}

plot <- ggplot(melt(model[ ,"plaq_density", ]), aes(x = Var2, y = value, group = factor(Var1))) +
  geom_line(aes(color = Var1), size = 1) +
  labs(x = "Cycle", y = "Plaque Density", color = "Individual") +
  geom_hline(yintercept = 0.82, linetype = "dashed")
plot


return(list(model_output = model, model_plot = plot))

}
