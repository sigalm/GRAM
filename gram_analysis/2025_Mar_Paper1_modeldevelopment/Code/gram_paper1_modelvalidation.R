# ********************************************************************************
# ======= GRAM PAPER 1: GRAM MODEL VALIDATION AND CALIBRATION ======= 
# ********************************************************************************

# TECHNICAL ##### 
# run only once, do not delete
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
source("gram_calibration/gram_benchmarking_helpers.r")
library(tableone)
sample1 <- readRDS("gram_data/acs_data/acs_sample_1.rds")
# sample2 <- readRDS("gram_data/acs_data/acs_sample_2.rds")
# sample3 <- readRDS("gram_data/acs_data/acs_sample_3.rds")
# ***************


# CALIBRATION ####
l.inputs1 <- l.inputs
l.inputs1[["n.cycle"]] <- 51
# l.inputs1[["m.lifetable"]] <- l.inputs[["m.lifetable"]] * 0.8
# 
# l.inputs1[["hr.mort_mci"]] <- 2
# l.inputs1[["hr.mort_mil"]] <- 3.5
# l.inputs1[["hr.mort_mod"]] <- 4.1
# l.inputs1[["hr.mort_sev"]] <- 10

# l.inputs1[["hr.mort_mil_age"]] <- l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(2, 2, 2) 
l.inputs1[["r.CDRslow_mean"]] <-  (seq(0, 1, length.out = 51)^1.5) * (1.5 *l.inputs[["r.CDRslow_mean"]])
# l.inputs1[["r.CDRfast_mean"]] <-  (seq(0.5, 1.5, length.out = 51)^2) * (2 * l.inputs[["r.CDRfast_mean"]])
l.inputs1[["p.HCARE_start"]] <- c(0.25,0.75)
l.inputs1[["n.ind"]] <- 50000
l.inputs1[["hr.mort_mci_age"]] <- c(1, 1, 1)
# l.inputs1[["hr.mort_mil_age"]] <- l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(2, 2, 2) 
l.inputs1[["hr.mort_mod_age"]] <- l.inputs1[["hr.mort_sev_age"]] <- c(1, 1, 1) 
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]] * 2.2 /4

mean(l.inputs1[["r.CDRslow_mean"]])
mean(l.inputs1[["r.CDRfast_mean"]])

sim_calib <- run_benchmarking(l.inputs1, "", sample1)

# Tables & Figures ####
## Table 1: Model Parameters ####
# created externally


## Table 2: Cohort Characteristics ####
tableone_vars <- c("AGE","SEX","RACEETH","EDU","INCOME","MEDBUR","HCARE","APOE4","SYN")
factor_vars <- list(
  SEX  = c("1" = "Male", "2" = "Female"),
  RACEETH   = c("0" = "Non-Hispanic White", "1" = "Non-Hispanic Black", "2" = "Hispanic"),
  INCOME = c("0" = "Low Income", "1" = "Medium Income", "2" = "High Income"),
  HCARE = c("0" = "No regular healthcare provider", "1" = "Has regular healthcare provider"),
  SYN  = c("0" = "Cognitively Healthy", "1" = "Cognitively Impaired"),
  APOE4  = c("0" = "Not Carrier", "1" = "Carrier")
)

cycle1 <- as.data.frame(t(sim_calib$sim$output[1, , ]))  %>% # transpose output so rows are indivs and cols are attrs
  mutate(across(names(factor_vars), ~ factor(.x, levels = names(factor_vars[[cur_column()]]), 
                                             labels = factor_vars[[cur_column()]])))

table2 <- CreateTableOne(vars = tableone_vars, data = cycle1, factorVars = names(factor_vars)) %>%
  print(showAllLevels = TRUE, printToggle = FALSE)

## Figure 1: Natural history of cognitive impairment ####
figure1 <- sim_calib$sim$fig.progression
figure1 <- figure1 + labs(title = NULL, subtitle = NULL)
ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figure1.pdf", plot = figure1, width = 10, height = 6, dpi = 300)

## Figure 2: Prevalence Results ####
figure2 <- compare_prevalence(sim_calib$sim, description = "", n = l.inputs1["n.ind"])
figure2 <- figure2$fig_prev_by_age + labs(title = NULL, subtitle = NULL)
ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figure2.pdf", plot = figure2, width = 10, height = 6, dpi = 300)

## Figure 3: Mortality Results ####
figure3 <- compare_mortality(sim_calib$sim, description = "", n = l.inputs1["n.ind"])
figure3 <- figure3$cumulative_mort[[1]] + labs(title = NULL, subtitle = NULL)
ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figure3.pdf", plot = figure3, width = 10, height = 6, dpi = 300)

## Figure 4: Reside Time Results ####
figure4 <- compare_reside_time(sim_calib$sim, description = "", n = l.inputs1["n.ind"])
figure4 <- figure4$fig_reside_time + labs(title = NULL, subtitle = NULL)
ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figure4.pdf", plot = figure4, width = 10, height = 6, dpi = 300)

# Reporting Results ####
# instructions: create variables to hold the following results that I will report in the main text. No need to make fancy figures/tables
# this is just so I have an easy way to see the number I want to report.
# 1. proportion with any cognitive impairment (i.e., MCI + dementia) at age 75
prop_CI_75 <- sum(sim_calib$sim$aggregated_results_totpop$state_trace[26, c("mci", "mil", "mod", "sev")])
prop_CI_75

# 2. proportion of healthy, MCI, and dementia at age 85
prop_CI_85 <- sim_calib$sim$aggregated_results_totpop$state_trace[36, ]
prop_CI_85

# 3. prevalence of MCI and dementia at age 85 (i.e., proportion among alive)
prev_CI_85 <- sim_calib$prev$dat %>%
  filter(age == 85)
prev_CI_85

# 4. Mean age of onset
sim_calib$age_onset$dat

# 5. avg absolute difference for prevalence between model vs benchmark at benchmark ages
benchmark_prevalence <- benchmark_prev_by_age %>%
  filter(condition == "dem", source == "Manly") %>%
  rename(prev_benchmark = prev)

deviation_prevalence <- sim_calib$prev$dat %>%
  filter(condition == "dem") %>%
  rename(prev_model = prev) %>%
  inner_join(benchmark_prevalence, by = c("age")) %>%
  select(age, raceeth, prev_model, prev_benchmark) %>%
  mutate(difference_abs = prev_model - prev_benchmark,
         difference_pct = (prev_model - prev_benchmark) / prev_benchmark)

avg_abs_difference_prevalence <- mean(abs(deviation_prevalence$difference_abs))
avg_abs_difference_prevalence

# 6. predicted annual mortality (<75 and 75+)
deviation_mortality <- sim_calib$mort$dat %>%
  select(age, model_rate_1000, benchmark_rate_1000, residual, cum_model_rate, cum_benchmark_rate, cum_residual) %>%
  mutate(residual_pct = residual / benchmark_rate_1000,
         cum_residual_pct = cum_residual / cum_benchmark_rate) %>%
  filter(is.finite(residual))

avg_abs_difference_annual_mortality_75under <- mean(abs(deviation_mortality$residual[deviation_mortality$age < 75])) 
avg_abs_difference_annual_mortality_75plus <- mean(abs(deviation_mortality$residual[deviation_mortality$age >= 75])) 

avg_abs_difference_annual_mortality_75under
avg_abs_difference_annual_mortality_75plus

# 7. predicted cumulative mortality (age 65, 75, 85)
cum_mort_by_age <- deviation_mortality %>%
  filter(age %in% c(65,75,85)) %>%
  select(age, cum_model_rate, cum_benchmark_rate, residual)
cum_mort_by_age

# 8. avg absolute difference for cumulative mortality between model vs benchmark
avg_abs_difference_cum_mortality <- mean(abs(deviation_mortality$cum_residual))
avg_abs_difference_cum_mortality

# 9. max absolute difference for cumulative mortality between model vs benchmark
max_abs_difference_cum_mortality <- deviation_mortality %>%
  filter(cum_residual == max(cum_residual))
max_abs_difference_cum_mortality 

# 9. median (IQR) reside time in MCI among those who ever develop MCI
mci_reside_times <- colSums(sim_calib$sim$output[,"SEV",] == 0, na.rm = TRUE)
median_time_mci <- summary(mci_reside_times[mci_reside_times > 0])
median_time_mci

# 10. median (IQR) reside time in dementia among those who ever develop dementia
dem_reside_times <- colSums(sim_calib$sim$output[,"SEV",] > 0, na.rm = TRUE)
median_time_dem <- summary(dem_reside_times[dem_reside_times > 0])
median_time_dem

# 11. overall life expectancy at age 50
sim_calib$sim$aggregated_results_totpop$alive.sum

# SUPPLEMENT ####
# Figure S1: Prevalence by risk factors

figS1_raceeth <- stratify_prevalence_by(sim_calib$sim, "RACEETH", 
                                        strat_labels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic"))
figS1_raceeth$plot_prev


figS1_income <- stratify_prevalence_by(sim_calib$sim, "INCOME", strat_labels = c("Low", "Medium", "High"))
figS1_income$plot_prev

figS1_edu <- stratify_prevalence_by(sim_calib$sim, "EDU", 
                                    strat_labels = c("Less than HS", "HS + Some College", "College or More"),
                                    strat_cutoffs = c(0, 12, 16, Inf))
figS1_edu$plot_prev

# figS1_medbur <- stratify_prevalence_by(sim_calib$sim, "MEDBUR", 
#                                     strat_labels = c("None", "Single Comorbidity", "Multiple Comorbidity"),
#                                     cutoffs = c(0, 1, 2, Inf))
# figS1_medbur$plot_prev

figS1_apoe4 <- stratify_prevalence_by(sim_calib$sim, "APOE4", c("Not Carrier", "Carrier"))
figS1_apoe4$plot_prev

# Combine all prevalence plots into one figure (e.g., 2 rows x 2 columns)
figS1_all <- 
  figS1_raceeth$plot_prev + labs(title = "(A) Race/Ethnicity") +
  figS1_income$plot_prev + labs(title = "(B) Income Group") +
  figS1_edu$plot_prev + labs(title = "(C) Educational Attainment") +
  figS1_apoe4$plot_prev + labs(title = "(D) APOE4 Carrier Status") +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

figS1_all

ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figureS1_combined.pdf", plot = figS1_all, width = 12, height = 14, dpi = 300)


# Figure S2: Age of onset by risk factors

# Function to compute per-group histogram proportions (each facet sums to 1)
make_hist_df <- function(df, group_var, binwidth = 3) {
  group_var <- rlang::ensym(group_var)
  hist_df <- df %>%
    filter(!is.na(age_onset), !is.na(!!group_var)) %>%
    group_by(!!group_var) %>%
    mutate(bin = cut(age_onset, breaks = seq(floor(min(age_onset)), ceiling(max(age_onset)), by = binwidth), right = FALSE)) %>%
    group_by(!!group_var, bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(!!group_var) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()
  # Add a bin_label column with more reader-friendly labels
  if (!is.null(hist_df$bin) && !all(is.na(hist_df$bin))) {
    bin_strs <- as.character(hist_df$bin)
    # Extract numeric bounds from bin labels (escape brackets and parens correctly)
    bin_label <- gsub("[\\[\\(]", "", bin_strs) # Remove [ or (
    bin_label <- gsub("\\]", "", bin_label) # Remove ]
    bin_label <- gsub("\\)", "", bin_label) # Remove )
    bin_label <- trimws(bin_label)
    # Split by comma and format as 'a-b'
    bin_label <- sapply(strsplit(bin_label, ","), function(x) paste0(as.integer(x[1]), "-", as.integer(as.numeric(x[2])-1)))
    hist_df$bin_label <- factor(bin_label, levels=unique(bin_label))
  } else {
    hist_df$bin_label <- NA
    message("[make_hist_df] Warning: bin_label was not created because 'bin' is NULL or all NA.")
  }
  print("[make_hist_df] Columns in output: ")
  print(colnames(hist_df))
  hist_df
}

onset_df <- data.frame(
  age_onset = sim_calib$sim$aggregated_results_totpop$age_at_onset,
  RACEETH = factor(sim_calib$sim$output[1, "RACEETH", ], levels = c(0, 1, 2), labels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic")),
  INCOME = factor(sim_calib$sim$output[1, "INCOME", ], levels = c(0, 1, 2), labels = c("Low Income", "Medium Income", "High Income")),
  EDU_BIN = cut(sim_calib$sim$output[1,"EDU",], breaks = c(0,12,16,Inf), labels = c("Less than HS", "HS + Some College", "College or More"), include.lowest = TRUE, right = FALSE),
  APOE4 = factor(sim_calib$sim$output[1, "APOE4", ], levels = c(0, 1), labels = c("Not Carrier", "Carrier"))
)

# S2 plots using precomputed proportions (each facet sums to 1)
hist_raceeth <- make_hist_df(onset_df, RACEETH)
hist_raceeth <- hist_raceeth[!is.na(hist_raceeth$bin_label) & hist_raceeth$bin_label != "NA-NA", ]
hist_raceeth$bin_label <- droplevels(factor(hist_raceeth$bin_label))
hist_income <- make_hist_df(onset_df, INCOME)
hist_income <- hist_income[!is.na(hist_income$bin_label) & hist_income$bin_label != "NA-NA", ]
hist_income$bin_label <- droplevels(factor(hist_income$bin_label))
hist_edu <- make_hist_df(onset_df, EDU_BIN)
hist_edu <- hist_edu[!is.na(hist_edu$bin_label) & hist_edu$bin_label != "NA-NA", ]
hist_edu$bin_label <- droplevels(factor(hist_edu$bin_label))
hist_apoe4 <- make_hist_df(onset_df, APOE4)
hist_apoe4 <- hist_apoe4[!is.na(hist_apoe4$bin_label) & hist_apoe4$bin_label != "NA-NA", ]
hist_apoe4$bin_label <- droplevels(factor(hist_apoe4$bin_label))

figS2_raceeth <- ggplot(hist_raceeth, aes(x = bin_label, y = prop)) +
  geom_col(fill = 'steelblue', color = 'white') +
  facet_wrap(~RACEETH, scales = "free_x") +
  labs(title = 'Age of Onset by Race/Ethnicity', x = 'Age of Onset', y = 'Proportion') +
  scale_x_discrete(drop = TRUE) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figS2_income <- ggplot(hist_income, aes(x = bin_label, y = prop)) +
  geom_col(fill = 'steelblue', color = 'white') +
  facet_wrap(~INCOME, scales = "free_x", drop = FALSE) +
  labs(title = 'Age of Onset by Income', x = 'Age of Onset', y = 'Proportion') +
  scale_x_discrete(drop = TRUE) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figS2_edu <- ggplot(hist_edu, aes(x = bin_label, y = prop)) +
  geom_col(fill = 'steelblue', color = 'white') +
  facet_wrap(~EDU_BIN, scales = "free_x") +
  labs(title = 'Age of Onset by Education', x = 'Age of Onset', y = 'Proportion') +
  scale_x_discrete(drop = TRUE) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

figS2_apoe4 <- ggplot(hist_apoe4, aes(x = bin_label, y = prop)) +
  geom_col(fill = 'steelblue', color = 'white') +
  facet_wrap(~APOE4, scales = "free_x") +
  labs(title = 'Age of Onset by APOE4 Carrier Status', x = 'Age of Onset', y = 'Proportion') +
  scale_x_discrete(drop = TRUE) +
  theme_minimal(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


figS2_all <- figS2_raceeth + labs(title = "(A) Race/Ethnicity") +
  figS2_income + labs(title = "(B) Income Group") +
  figS2_edu + labs(title = "(C) Educational Attainment") +
  figS2_apoe4 + labs(title = "(D) APOE4 Carrier Status") +
  plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())

figS2_all

figS2_raceeth$data %>%
  group_by(RACEETH) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  select(RACEETH, bin, bin_label, count, prop)

figS2_income$data %>%
  group_by(INCOME) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  select(INCOME, bin, bin_label, count, prop)

figS2_edu$data %>%
  group_by(EDU_BIN) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  select(EDU_BIN, bin, bin_label, count, prop)

figS2_apoe4$data %>%
  group_by(APOE4) %>%
  slice_max(prop, n = 1, with_ties = FALSE) %>%
  select(APOE4, bin, bin_label, count, prop)
  
ggsave("gram_analysis/2025_Mar_Paper1_modeldevelopment/Figures/figureS2_combined.pdf", plot = figS2_all, width = 12, height = 14, dpi = 300)

