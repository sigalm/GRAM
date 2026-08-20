library(dplyr)

a.out <- sim_calib$sim$output


# Extract needed variables
n_cycles <- dim(a.out)[1]
n_ind <- dim(a.out)[3]

sex_vec <- a.out[1, "SEX", ]+1   # baseline sex: 1 = men, 2 = women

# Reshape simulation output to long format
age_long   <- as.vector(a.out[ , "AGE", ])
sex_long   <- rep(sex_vec, each = n_cycles)
alive_long <- as.vector(a.out[ , "ALIVE", ])
syn_long <- as.vector(a.out[ , "SYN", ])
sev_long   <- as.vector(a.out[ , "SEV", ])
memloss_long <- as.vector(a.out[ , "MEMLOSS", ])
dx_long <- as.vector(a.out[, "DX", ])

# Define dementia: SEV == 1,2,3
dementia_long <- (sev_long %in% c(1, 2, 3)) & (alive_long == 1)
mci_long <- c(sev_long == 0) & (alive_long == 1) & !memloss_long
tci_long <- c(syn_long == 0.5) & (alive_long == 1)

# Data frame for analysis (includes age)
df <- data.frame(
  age = age_long,
  sex = factor(sex_long, levels = c(1, 2), labels = c("Men", "Women")),
  alive = alive_long,
  tci = tci_long,
  memloss = memloss_long,
  mci = mci_long,
  dementia = dementia_long,
  dx = dx_long
)

# Restrict to alive person-cycles
df_alive <- df[df$alive == 1, ]

# --- Prevalence for full population ---
cases_by_sex <- df_alive %>%
  filter(dementia == TRUE) %>%
  group_by(sex) %>%
  summarize(
    total_dementia = n(),
    .groups = "drop"
  ) %>%
  mutate(
    proportion = total_dementia / sum(total_dementia)
  )

print(cases_by_sex)

# --- Prevalence for age 65+ ---
prevalence_by_sex_65plus <- df_alive %>%
  filter(age >= 65) %>%
  group_by(sex) %>%
  summarize(
    total_dementia = sum(dementia),
    total_alive = n(),
    prevalence = total_dementia / total_alive
  )

print(prevalence_by_sex_65plus)


(prevalence_65 <- df_alive %>%
  filter(age == 65, dx == 0) %>%
  summarize(tci = sum(tci)/n(),
            memloss = sum(memloss)/n(),
            mci = sum(mci, na.rm = TRUE)/n(),
            dem = sum(dementia)/n()))


(prevalence_70 <- df_alive %>%
  filter(age == 70, dx == 0) %>%
    summarize(tci = sum(tci)/n(),
              memloss = sum(memloss)/n(),
              mci = sum(mci, na.rm = TRUE)/n(),
              dem = sum(dementia)/n()))

(prevalence_80 <- df_alive %>%
    filter(age == 80, dx == 0) %>%
    summarize(tci = sum(tci)/n(),
              memloss = sum(memloss)/n(),
              mci = sum(mci, na.rm = TRUE)/n(),
              dem = sum(dementia)/n()))


