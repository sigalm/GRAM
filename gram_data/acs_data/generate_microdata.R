## Generate microdata

install.packages("ipumsr")
library(ipumsr)
library(dplyr)
set.seed(20240324)

# Load the IPUMS DDI metadata
ddi <- read_ipums_ddi("gram_data/acs_data/usa_00004.xml")

# Read the microdata file
acs_data_raw <- read_ipums_micro(ddi)

# Select relevant variables
acs_data <- acs_data_raw %>%
  select(AGE, SEX, RACE, HISPAN, HCOVANY, EDUCD, FTOTINC, DIFFREM, PERWT)

# Clean data
acs_data <- acs_data %>%
  rename(
    EDUC = EDUCD,
    INCOME = FTOTINC,
    RACE = RACE,
    HISP = HISPAN,
    INSURANCE = HCOVANY,
    COGCON = DIFFREM
  ) %>%
  mutate(
    EDUC = case_when(
      EDUC %in% c(000, 001, 002) ~ 0,
      EDUC == 010 ~ 2.5,  # Nursery-Grade 4 (Unspecified)
      EDUC == 011 ~ 0,    # Nursery, preschool
      EDUC == 012 ~ 0,    # Kindergarten
      EDUC == 013 ~ 2.5,  # Grade 1-4 (Unspecified)
      EDUC == 014 ~ 1,
      EDUC == 015 ~ 2,
      EDUC == 016 ~ 3,
      EDUC == 017 ~ 4,
      EDUC == 020 ~ 6.5,  # Grades 5-8 (Unspecified)
      EDUC == 021 ~ 5.5,  # Grades 5-6 (Unspecified)
      EDUC == 022 ~ 5,
      EDUC == 023 ~ 6,
      EDUC == 024 ~ 7.5,  # Grades 7-8 (Unspecified)
      EDUC == 025 ~ 7,
      EDUC == 026 ~ 8,
      EDUC == 030 ~ 9,
      EDUC == 040 ~ 10,
      EDUC == 050 ~ 11,
      EDUC %in% c(060, 061, 062, 063, 064) ~ 12,  # HS Grad/GED
      EDUC == 065 ~ 13,  # Some college, no degree
      EDUC %in% c(070, 071) ~ 14,  # 1+ years college
      EDUC %in% c(080, 081, 082, 083) ~ 14,  # Associate’s degree
      EDUC == 090 ~ 15,
      EDUC %in% c(100, 101) ~ 16,  # Bachelor's degree
      EDUC %in% c(110, 111, 112, 113) ~ 17,  # 5-8 Years College
      EDUC == 114 ~ 18,  # Master’s
      EDUC == 115 ~ 19,  # Professional
      EDUC == 116 ~ 20,  # Doctorate
      EDUC == 999 ~ NA  # Missing values
    ),
    HISP = case_when(
      HISP == 0 ~ 0,  # Not Hispanic
      HISP %in% c(1,2,3,4) ~ 1,  # Mexican, Puerto Rican, Cuban, Other
      HISP == 9 ~ NA  # Missing
    ),
    RACEETH = case_when(
      RACE == 1 & HISP == 0 ~ 0,  # Non-Hispanic White 
      RACE == 2 & HISP == 0 ~ 1,  # Non-Hispanic Black
      HISP == 1 ~ 2 # Hispanic (can be other race too)
    ),
    INCOME_CAT = case_when(
      INCOME < 0 ~ NA,
      INCOME >= 9999998 ~ NA,
      INCOME < 9000 ~ 0,
      INCOME < 36000 ~ 1,
      INCOME >= 36000 ~ 2
    ),
    INSURANCE = INSURANCE - 1) %>%
  filter(AGE == 50, !is.na(RACEETH))

head(acs_data)
summary(acs_data$EDUC)
table(acs_data$EDUC, useNA = "ifany")
table(acs_data$INCOME_CAT, useNA = "always")

## Add medical burden
# See https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-023-15768-8 figure 1
prev.no_of_conditions_male <- c(0.155, 0.205, 0.220, 0.175, 0.105, 0.060, 0.045, 0.020, 0.018, 0.005, 0.001)
prev.no_of_conditions_female <- c(0.175, 0.205, 0.190, 0.175, 0.140, 0.075, 0.040, 0.018, 0.015, 0.002, 0.002)

prev.no_of_conditions <- (prev.no_of_conditions_female + prev.no_of_conditions_male) / 2

set.seed(123)  

oddsratio.2plus_conditions <- case_when(
  acs_data$EDUC >= 16 ~ 1,
  acs_data$EDUC >= 12 ~ 1.32,
  acs_data$EDUC < 12  ~ 1.58) 

prev_2plus <- sum(prev.no_of_conditions) - sum(prev.no_of_conditions[1:2])  
prev_2plus_col <- prev_2plus / (prop.table(table(oddsratio.2plus_conditions))[1] + 
                                  1.32 * prop.table(table(oddsratio.2plus_conditions))[2] +
                                  1.58 * prop.table(table(oddsratio.2plus_conditions))[3])

# Prevalence of 2+ conditions among college grads is 0.534

p.2plus_conditions <- oddsratio.2plus_conditions * prev_2plus_col

# Generate binary indicator for having 2+ conditions
tmp_rand <- runif(n = nrow(acs_data))
tmp_2plus_conditions <- qbinom(p = tmp_rand, size = 1, prob = p.2plus_conditions)

acs_data$MEDBUR <- case_when(
  tmp_2plus_conditions == 1 ~ sample(2:10, nrow(acs_data), replace = TRUE, 
                                     prob = prev.no_of_conditions[3:11] / sum(prev.no_of_conditions[3:11])),
  tmp_2plus_conditions == 0 ~ sample(0:1, nrow(acs_data), replace = TRUE, 
                                     prob = prev.no_of_conditions[1:2] / sum(prev.no_of_conditions[1:2]))
)


## Add APOE4 data (random, 25% prevalence)
acs_data$APOE4 <- rbinom(n = nrow(acs_data), size = 1, prob = 0.25)

# Check and save full data
table(acs_data$RACEETH, acs_data$HISP)
table(acs_data$AGE)
prop.table(table(acs_data$SEX))

summary(acs_data) # Income has NAs, remove them

acs_data <- acs_data %>%
  filter(!is.na(INCOME_CAT))

saveRDS(acs_data, "gram_data/acs_data/acs_age50.RDS")


# Create sample cohorts: 3 sets with 10,000 people, accounting for person weights.
generate_synthetic_sample <- function(pop_data, target_size, seed) {
  set.seed(seed)  # Ensure reproducibility for each dataset
  pop_data %>%
    slice_sample(n = target_size, weight_by = PERWT, replace = TRUE)
}

sample1 <- generate_synthetic_sample(acs_data, target_size = 10000, seed = 1001)
sample2 <- generate_synthetic_sample(acs_data, target_size = 10000, seed = 1002)
sample3 <- generate_synthetic_sample(acs_data, target_size = 10000, seed = 1003)


# Validation checks
colSums(is.na(acs_data))  

sapply(list(sample1, sample2, sample3), nrow)

summary(sample1$EDUC)
summary(sample2$EDUC)
summary(sample3$EDUC)

prop.table(table(sample1$SEX))
prop.table(table(sample2$SEX))
prop.table(table(sample3$SEX))

prop.table(table(sample1$APOE4))
prop.table(table(sample2$APOE4))
prop.table(table(sample3$APOE4))

summary(acs_data$AGE)
summary(acs_data$EDUC)
summary(acs_data$MEDBUR)
summary(acs_data$APOE4)
prop.table(table(acs_data$APOE4))
prop.table(table(acs_data$MEDBUR))


# Save samples
saveRDS(sample1, "gram_data/acs_data/acs_sample_1.rds")
saveRDS(sample2, "gram_data/acs_data/acs_sample_2.rds")
saveRDS(sample3, "gram_data/acs_data/acs_sample_3.rds")


