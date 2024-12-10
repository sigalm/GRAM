# *****************************************************
# ======= CALCULATE DEMENTIA-ADJUSTED MORTALITY ======= 
# *****************************************************

library(tidyverse)

path_allcause <- "../gram_data/allcause_mortality_by_age.txt"
path_dementia <- "../gram_data/dementia_mortality_by_age.txt"

allcause_raw <- read_tsv(path_allcause, n_max = 101)  # Crude rates of death per 100,000 by age
dementia_raw <- read_tsv(path_dementia, n_max = 60)   # Crude rates of death per 1,000 by age

str(allcause_raw)

allcause <- allcause_raw %>%
  janitor::clean_names() %>%
  select(-notes, -single_year_ages) %>%
  rename(age = single_year_ages_code,
         rate = crude_rate) %>%
  mutate(age = as.numeric(age),
         population = as.numeric(population),
         rate = deaths/population) %>%       
  subset(age >= 50)



str(allcause)      # all converted to numeric
summary(allcause)  # no NAs in age, 16 in rate (those aged 85+), all rates are non-zero



str(dementia_raw)

dementia <- dementia_raw %>%
  select(-Notes, -`Single-Year Ages`) %>%
  rename(age = `Single-Year Ages Code`,
         deaths = Deaths,
         population = Population,
         rate = `Crude Rate`) %>%
  mutate(age = as.numeric(age),
         population = as.numeric(population),
         rate = deaths/population) %>%        
  subset(age >= 50)



str(dementia)      # all converted to numeric
summary(dementia)  # no NAs in age, 16 in rate (those age 85+), all rates are non-zero



# rates for very young and very old ages are unavailable, extrapolate using available data

plot(allcause$age, allcause$rate)  # visualize rate by age -- it is exponential
plot(dementia$age, dementia$rate)  # visualize rate by age -- it is exponential

# Since rate grows exponentially, will use log-transformed linear model

any(allcause$rate==0, na.rm=TRUE) # check to see if any zeros (log(0) is problematic) -- there are none
any(dementia$rate==0, na.rm=TRUE) # check to see if any zeros (log(0) is problematic) -- there are none


# Extrapolate all-cause mortality 
summary(allcause_log <- lm(log(rate) ~ age, data = na.omit(allcause[allcause$age>=76, ])))   # using all data points gave weird fits, so only using ages 76 and up to predict
# good fit (R^2 = 0.99), p(age) ~ 0

allcause_predict <- data.frame(age = allcause$age[is.na(allcause$rate)])
allcause_predict$rate_predicted <- exp(predict(allcause_log, newdata = allcause_predict))

allcause_clean <- allcause %>%
  left_join(allcause_predict, by = "age") %>% 
  mutate(rate = coalesce(rate, rate_predicted)) %>%
  select(-rate_predicted)
summary(allcause_clean)  # no NAs in rate2 and probability2

plot(allcause_clean$age, allcause_clean$rate) # smooth exponential curve

# DISCARD 
# # Extrapolate dementia mortality 
# summary(dementia_log <- lm(log(rate) ~ age, data = na.omit(dementia)))   # good fit (R^2 ~= 1), p(age) ~ 0
# 
# dementia_predict <- data.frame(age = dementia$age[is.na(dementia$rate)])
# dementia_predict$rate_predicted <- exp(predict(dementia_log, newdata = dementia_predict))  # good fir
# 
# dementia_clean <- dementia %>%
#   left_join(dementia_predict, by = "age") %>% 
#   mutate(rate = coalesce(rate, rate_predicted)) %>%
#   select(-rate_predicted)
# 
# summary(dementia_clean)
# plot(dementia_clean$age, dementia_clean$rate) # smooth exponential curve
# 
# 
# # Now that we have all-cause and dementia-specific mortality rates by age, calculate non-dementia mortality
# non_dementia <- allcause_clean[ ,c("age","rate")] %>%
#   right_join(dementia_clean[ ,c("age","rate")], by="age") %>%
#   rename(rate_allcause = rate.x,
#          rate_dementia = rate.y) %>%
#   mutate(rate_non_dementia = pmax(rate_allcause - rate_dementia, 1e-10))
# 
# summary(non_dementia)
# 
# plot(non_dementia$age, non_dementia$rate_allcause)
# plot(non_dementia$age, non_dementia$rate_dementia)  
# plot(non_dementia$age, non_dementia$rate_non_dementia)  
# 


# Or extrapolate after subtracting -- might be better
non_dementia2 <- allcause_clean[ ,c("age","rate")] %>%
  right_join(dementia[ ,c("age","rate")], by="age") %>%
  rename(rate_allcause = rate.x,
         rate_dementia = rate.y) %>%
  mutate(rate_non_dementia = pmax(rate_allcause - rate_dementia, 1e-10))

summary(non_dementia2)

summary(non_dementia_log <- lm(log(rate_non_dementia) ~ age, data = na.omit(non_dementia2)))   # good fit (R^2 ~= 1), p(age) ~ 0

non_dementia2_predict <- data.frame(age = non_dementia2$age[is.na(non_dementia2$rate_non_dementia)])
non_dementia2_predict$rate_predicted <- exp(predict(non_dementia_log, newdata = non_dementia2_predict))  # good fir

non_dementia2_clean <- non_dementia2 %>%
  left_join(non_dementia2_predict, by = "age") %>% 
  mutate(rate_non_dementia = coalesce(rate_non_dementia, rate_predicted)) %>%
  select(-rate_predicted, -rate_dementia)

summary(non_dementia2_clean)
plot(non_dementia2_clean$age, non_dementia2_clean$rate_non_dementia) # smooth exponential curve

non_dementia2_clean <- non_dementia2_clean %>%
  mutate(prob_allcause = 1-exp(-rate_allcause),
         prob_non_dementia = 1-exp(-rate_non_dementia))

non_dementia2_grouped <- non_dementia2_clean %>%
  mutate(age_group = cut(age, breaks = seq(50, 100, by = 5), right = FALSE)) %>%
  group_by(age_group) %>%
  summarise(avg_rate_allcause = mean(rate_allcause),
            avg_rate_non_dementia = mean(rate_non_dementia)) %>%
  mutate(avg_prob_allcause = 1-exp(-avg_rate_allcause),
         avg_prob_non_dementia = 1-exp(-avg_rate_non_dementia))



# # Convert rates to probabilities
# non_dementia_probs <- non_dementia %>%
#   mutate(p_death_non_dementia = 1 - exp(-rate_non_dementia)) %>%
#   select(-rate_allcause,-rate_dementia,-rate_non_dementia)

saveRDS(non_dementia2_clean, "../gram_data/non_dementia_mortality_prob_by_age_v2.RDS")


# 
# # Alternative method, using Jim's spreadsheet
# non_dementia2_grouped <- non_dementia2_clean %>%
#   mutate(age_group = cut(age, breaks = seq(50, 100, by = 5), right = FALSE)) %>%
#   group_by(age_group) %>%
#   summarise(avg_rate_allcause = mean(rate_allcause),
#             avg_rate_dementia = mean(rate_dementia),
#             avg_rate_non_dementia = mean(rate_non_dementia)) %>%
#   mutate(avg_prob_allcause = 1-exp(-avg_rate_allcause),
#          avg_prob_dementia = 1-exp(-avg_rate_dementia),
#          avg_prob_non_dementia = 1-exp(-avg_rate_non_dementia))


prevalence_grouped <- as.data.frame(scenario1.1$aggregated_results_totpop$state_trace) %>%
  mutate(age = 50:99) %>%
  mutate(alive = 1 - dth,
         age_group = cut(age, breaks = seq(50, 100, by = 5), right = FALSE)) %>%
  mutate(healthy = healthy / alive,
         mci = mci / alive,
         mil = mil / alive,
         mod = mod / alive,
         sev = sev / alive) %>%
  group_by(age_group) %>%
  summarize(avg_healthy = mean(healthy),
            avg_mci = mean(mci),
            avg_mil = mean(mil),
            avg_mod = mean(mod),
            avg_sev = mean(sev))


lifetable <- read.csv("../gram_data/table_prob_die_next_year.csv")
lifetable_grouped <- lifetable %>%
  mutate(age_group = cut(age, breaks = seq(50, 100, by = 5), right = FALSE)) %>%
  group_by(age_group) %>%
  summarise(avg_qx = mean(qx))
