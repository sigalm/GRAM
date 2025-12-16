# Predict baseline incidence rates (per 1000 person-years)
# Based on Gillis et al 2019: https://pmc.ncbi.nlm.nih.gov/articles/PMC6416157, Table 1


# Use mid-points of age groups
df <- data.frame(
  age = c(67, 72,   77,   82,   87),
  mci = c(14, NA, 22.5, 40.9, 60.1)  
)

# Relationship is exponential. Will do a log-transformed linear model
summary(incidence_log <- lm(log(mci) ~ age, data = na.omit(df)))

# Good fit, p<0.05 for age coefficient.

# Predict incidence rate for ages 50 to 100
df_predict <- data.frame(
  age = 50:100,
  mci = rep(NA, 51))

mci_predicted <- exp(predict(incidence_log, newdata = df_predict)) 

# Put all in new dataframe to feed into model as baseline hazard and save 
df_predict <- df_predict %>%
  mutate(mci = mci_predicted)

plot(df_predict$mci)

saveRDS(df_predict, "gram_data/mci_incidence/mci_incidence_rate_by_age.RDS")
