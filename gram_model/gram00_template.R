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
l.inputs1[["scenario"]] <- "Natural history, no intervention"  # always start by changing the scenario title
l.inputs1[["r.CDRslow_mean"]] <- 0.3
# ...

# Run model with updated inputs
scenario1 <- f.wrap_run(l.inputs1)

# Save and format outputs you need as separate variables for easy access
fig_scenario1 <- scenario1$fig.progression$fig.progression_true
reside_time_scenario1 <- format_reside_time_table(scenario1$aggregated_results_totpop$reside_time_by_age_of_onset)
prevalence_scenario1 <- format_prevalence_table(scenario1$aggregated_results_totpop$prevalence_grouped)

## SCENARIO 2 NAME ####
# repeat the section above for as many scenarios as needed


# EXPORT RESULTS ####

# Export figures (change file path extension if you want a different file type e.g., pdf, jpeg)
ggsave("[replace with file path].png", plot = fig_scenario1, width = 10, height = 6, dpi = 300)

## Export tables (you can add multiple tables to a single word doc by including multiple tables in the function call
save_as_docx("[replace with table title]" = reside_time_scenario1, path = "[replace with file path]")

