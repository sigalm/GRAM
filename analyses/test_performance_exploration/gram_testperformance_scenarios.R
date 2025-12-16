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
## SCENARIO 1a: Strict test (high specificity, low sensitivity), cognitive concerns not considered ####

l.inputs1a <- l.inputs
scenario_strict_noconcerns <- f.wrap_run(l.inputs1a)

## SCENARIO 1b: Lax test (low specificity, high sensitivity), cognitive concerns not considered ####
l.inputs1b <- l.inputs
l.inputs1b[["sens_BHA"]] <- l.inputs[["m.BHA_performance"]]$sens[[2]]
l.inputs1b[["spec_BHA"]] <- l.inputs[["m.BHA_performance"]]$spec[2]

scenario_lax_noconcerns <- f.wrap_run(l.inputs1b)

## SCENARIO 2a: Strict test (high specificity, low sensitivity), spontaneous concerns ####

l.inputs2a <- l.inputs1a
l.inputs2a[["m.cogcon"]] <- l.inputs[["m.cogcon_spon"]]

scenario_strict_sponconcerns <- f.wrap_run(l.inputs2a)

## SCENARIO 2b: Lax test (low specificity, high sensitivity), spontaneous concerns ####
l.inputs2b <- l.inputs1b
l.inputs2b[["m.cogcon"]] <- l.inputs[["m.cogcon_spon"]]

scenario_lax_sponconcerns <- f.wrap_run(l.inputs2b)

## SCENARIO 3a: Strict test (high specificity, low sensitivity), elicited concerns ####

l.inputs3a <- l.inputs2a
l.inputs3a[["m.cogcon"]] <- l.inputs[["m.cogcon_elic"]]

scenario_strict_elicconcerns <- f.wrap_run(l.inputs3a)

## SCENARIO 3b: Lax test (low specificity, high sensitivity), elicited concerns ####

l.inputs3b <- l.inputs2b
l.inputs3b[["m.cogcon"]] <- l.inputs[["m.cogcon_elic"]]

scenario_lax_elicconcerns <- f.wrap_run(l.inputs3b)


# REVIEW RESULTS ####
f.analyze_test_performance(scenario_strict_noconcerns, 25)
f.analyze_test_performance(scenario_strict_sponconcerns, 25)
f.analyze_test_performance(scenario_strict_elicconcerns, 25)

f.analyze_test_performance(scenario_lax_noconcerns, 25)
f.analyze_test_performance(scenario_lax_sponconcerns, 25)
f.analyze_test_performance(scenario_lax_elicconcerns, 25)

## both tests negative, 1 of the 2 tests positive, both tests positive (assuming independence)
## separate results by those who remained healthy, and those who remained in MCI
## repeat for strict and lax

## for the subset of people who are healthy at 70 and MCI at 72
## how many tested tn the tp, how many tested tn and fn, those tested fp and tp, those went back fp to fn

## in the Possin 2018 study people who got tested were people who had some concerns -- that's how they ended up in the study
## gold standard is a comprehensive evaluation

# do we need to consider if someone has a caregiver/family member to report a concern?
# do we test differently based on this?

## for someone who had an elicited concern at age 70, tested negative on the BHA, 
# how likely will that person to have either elicited or spontaneous concerns at age 85.
# this needs to be a function of their true disease state. 
# if they are cog normal, have a spont concern, but test negative, 2-5 years later they're less likely to have a concern.
# who really did have MCI and tested false neg
# doing the BHA result 

## incorporate access to healthcare also interacts with expressing concern -- you can't express concern if you don't go to doctor
## expressing concern is usually spontaneous, unprompted.
## normals are known to NOT have concerns

## what do we know about test performance among those with concerns -- not much

## Kate to send Kaiser data for real-world test performance using tabcat. - more diverse population

## Avg time between BHA tests is 1.5 years and there is data on the improvement 
## Look at Elena's study if updated test performance (she had larger sample size)

### OTHER NOTES ON GENERAL MODEL:
# Use the same random number for BHA over time e.g., always use a.random[1,"BHA", alive]


# # 
# 
# 
# # Run model with updated inputs
# scenario1 <- f.wrap_run(l.inputs1)
# 
# # Save and format outputs you need as separate variables for easy access
# fig_scenario1 <- scenario1$fig.progression$fig.progression_true
# reside_time_scenario1 <- f.format_reside_time_table(scenario1$aggregated_results_totpop$reside_time$noncensored)
# prevalence_scenario1 <- f.format_prevalence_table(scenario1$aggregated_results_totpop$prevalence_grouped)
# age_onset <- f.make_histogram(scenario1$aggregated_results_totpop$age_at_mci, lab = "Age at MCI Onset")
# mean(scenario1$aggregated_results_totpop$age_at_mci, na.rm = TRUE)







# EXPORT RESULTS ####

# Export figures (change file path extension if you want a different file type e.g., pdf, jpeg)
# ggsave("[replace with file path].png", plot = fig_scenario1, width = 10, height = 6, dpi = 300)

## Export tables (you can add multiple tables to a single word doc by including multiple tables in the function call
# save_as_docx("[replace with table title]" = reside_time_scenario1, path = "[replace with file path]")

