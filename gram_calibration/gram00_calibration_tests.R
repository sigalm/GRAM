# ********************************************************************************
# ======= GRAM SCENARIO RUN TEMPLATE [REPLACE WITH DESCRIPTIVE TITLE] ======= 
# ********************************************************************************

# TECHNICAL ##### 
# run only once, do not delete
source("gram_model/gram_01setup.r")
source("gram_model/gram_02helpers.r")
source("gram_model/gram_03simulation.r")
# ***************


# SCENARIOS ####
## SCENARIO 1 NAME (EXAMPLE: Natural progression, no treatment) ####

# Make a copy of inputs file to edit
l.inputs1 <- l.inputs

# Look up variable names in l.inputs
names(l.inputs)

# Update variables to define scenario
l.inputs1[["coef_MEDBUR"]] <- 0.2
l.inputs1[["amplification_MEDBUR"]] <- 0.025

l.inputs2 <- l.inputs
l.inputs2[["coef_MEDBUR"]] <- 0.1
l.inputs2[["amplification_MEDBUR"]] <- 0.025
l.inputs2[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]]*4


# values from sensitivity analysis model (as opposed to full model) from Table 3 of Angevaare et al. - this is worse!
l.inputs1[["log_EDU"]] <- log(0.96)
l.inputs1[["log_SEX"]] <- log(0.80)
l.inputs1[["log_RACEETHhisp"]] <- log(0.84)
l.inputs1[["log_APOE4"]] <- log(1.11)
l.inputs1[["log_MEDBUR"]] <- log(1.08)
l.inputs1[["log_INCOMEmed"]] <- log(0.83)
l.inputs1[["log_INCOMEhi"]] <- log(0.79)


# try reducing the underlying hazard rate (incidence rate) for MCI
l.inputs1[["m.hr_mci"]] <- l.inputs[["m.hr_mci"]]*4

# ...

# Run model with updated inputs
scenario1 <- f.wrap_run(l.inputs1)
l.inputs[["coef_MEDBUR"]]

scenario2 <- f.wrap_run(l.inputs2)
mean(scenario2$aggregated_results_totpop$age_at_mci, na.rm = TRUE)


# Save and format outputs you need as separate variables for easy access
fig_scenario1 <- scenario1$fig.progression$fig.progression_true
reside_time_scenario1 <- f.format_reside_time_table(scenario1$aggregated_results_totpop$reside_time$noncensored)
prevalence_scenario1 <- f.format_prevalence_table(scenario1$aggregated_results_totpop$prevalence_grouped)
age_onset <- f.make_histogram(scenario1$aggregated_results_totpop$age_at_mci, lab = "Age at MCI Onset")
mean(scenario1$aggregated_results_totpop$age_at_mci, na.rm = TRUE)
mean(scenario1$output[1,"MEDBUR",], na.rm = TRUE)
mean(scenario1$output[25,"MEDBUR",], na.rm = TRUE)
mean(scenario1$output[1,"EDU",], na.rm = TRUE)
hist(scenario1$output[1,"EDU",])
scenario1$output[,"MEDBUR",5]

## SCENARIO 2 NAME ####
# repeat the section above for as many scenarios as needed


# EXPORT RESULTS ####

# Export figures (change file path extension if you want a different file type e.g., pdf, jpeg)
ggsave("[replace with file path].png", plot = fig_scenario1, width = 10, height = 6, dpi = 300)

## Export tables (you can add multiple tables to a single word doc by including multiple tables in the function call
save_as_docx("[replace with table title]" = reside_time_scenario1, path = "[replace with file path]")

