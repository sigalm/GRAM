######################################## PAPER 2 WORD REPORT: TABLES + MAIN FIGURE ########################################
# Rebuilds every paper 2 table, and the counts-and-shares figure, from a run that has
# already happened. Nothing is re-simulated and nothing is re-derived here: the analysis
# script itself is sourced with `run_scenarios <- FALSE`, so the tables in the Word
# document are the same objects the analysis produces, not a second copy of the logic.
#
# Run from the project root. Set report_run_id to the timestamp of the run to report.

report_run_id <- "20261123_101738"

# Options, not variables: the analysis script sources model/setup.R, which clears the
# global environment before it reaches these.
options(gram.run_id        = report_run_id,
        gram.run_scenarios = FALSE,   # reuse that run's sim_results/ and test_perf_results/
        gram.make_figures  = FALSE,   # the jpegs already exist; this script writes its own PNG
        # Everything the paper and its supplement report. The analysis script's own
        # default is the main arms only, and a table whose arms are not loaded is skipped.
        gram.scenarios     = c("r1bhapos", "s1bhapos_emr", "u3bhapos_rand50",
                               "s1bhapos_question", "u3bhapos_nonrand",
                               "r1pcp_perfect", "s1pcp_perfect_emr", "u3pcp_perfect_rand50"))

source("analyses/testing_strategies/gram_paper2_testingstrategies.R")
report_run_id <- getOption("gram.run_id")

report_dir <- "analyses/testing_strategies/word_report"
fig_dir    <- file.path(report_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)


## Counts-and-shares figure, sized for the page rather than for a slide ####
# Text in a ggplot scales with the physical size it is RENDERED at, not the size it is
# displayed at, so the analysis script's 17 in wide figure, shown across a 9.2 in text
# block, would lose nearly half its type size. Rendering close to display size is what
# keeps the type near its set size on the page.
#
# The width is the page rather than the analysis script's: 15 in across the landscape
# text block. One row now, so 6.6 in high rather than the two rows' 9.6.
#
# Sizes are set for a figure that gets pasted into Word and shrunk: base_size drives the
# strip and axis titles, label_size the band labels, and the two are kept close so the
# labels are not what becomes unreadable first. right_pad is the gutter the labels sit
# in, a share of the x axis while the labels are sized in points, so it grows with them.
#
# The tested and test counts live in a table now, so no strategy_stats; the note the
# "Not eligible*" asterisk points at goes in the figure's caption in the document.
counts_and_shares <- plot_counts_and_shares(
  test_data_for(main_keys_by_intensity),
  scenario_names = labels_for(main_keys_by_intensity),
  ages            = plot_ages,
  base_size       = 16,
  label_size      = 4.8,
  band_gap        = 0.058,
  right_pad       = 0.8,     # room for "Not tested, impaired"
  # Untested split healthy / impaired under one teal fill, the impaired part hatched:
  # the misses the program never tested.
  hatch_missed    = TRUE,
  # The shares row alone: the counts row above it told the same story twice. The
  # Word heading is the figure title, so the plot carries none of its own.
  rows            = "shares",
  bottom_title    = NULL,
  bottom_subtitle = NULL,
  bottom_caption  = NULL)

ggsave(file.path(fig_dir, "counts_and_shares.png"), plot = counts_and_shares,
       width = 15, height = 6.6, dpi = 300, bg = "white")


## Annual performance figure, at the same page width as the counts-and-shares one ####
annual_performance <- plot_annual_performance(
  test_data_for(main_keys_by_intensity),
  scenario_names = labels_for(main_keys_by_intensity),
  ages = plot_ages, base_size = 20, label_size = 5.5)

ggsave(file.path(fig_dir, "annual_performance.png"), plot = annual_performance,
       width = 15, height = 6.6, dpi = 300, bg = "white")


## Test characteristics ####
# BHA_GS (model/test_properties.R) indexes sensitivity by the branch order in
# f.update_BHA(), not by severity: [1] TCI, [2] non-progressive memory loss, [3] MCI,
# [4] dementia (one value shared by all three dementia severities).
bha_params <- data.frame(
  Parameter = c("Sensitivity, MCI", "Sensitivity, dementia", "Specificity"),
  Value     = sprintf("%.2f", c(BHA_GS$sens[3], BHA_GS$sens[4], BHA_GS$spec)),
  check.names = FALSE)

bha_other_sens <- sprintf("%.2f", c(tci = BHA_GS$sens[1], memloss = BHA_GS$sens[2]))


## Table ordering ####
# Every table reads least to most intensive: Reactive, Selective, Inclusive, and for the
# PCP set too. The registry's own order is the figure order, which is not the same
# thing. Only the perfect PCP is reported (see the registry).
pcp_keys_ordered <- c("r1pcp_perfect", "s1pcp_perfect_emr", "u3pcp_perfect_rand50")

# A perfect PCP dismisses every positive given before onset, so the PCP arms detect no
# one before onset by construction, and a column of zeros says nothing.
drop_early <- function(x) x[!startsWith(names(x), "Coincidentally detected before onset")]


## Save ####
saveRDS(list(
  run_id            = report_run_id,
  n_ind             = l.inputs_calibrated[["n.ind"]],
  testing_window    = testing_window,
  end_age           = end_age,
  counts_and_shares = file.path("figures", "counts_and_shares.png"),
  uptake_inputs     = uptake_inputs_for(
    c("r1bhapos", "s1bhapos_emr", "s1bhapos_question", "u3bhapos_rand50", "u3bhapos_nonrand"),
    labels = c("Reactive", "Selective – eRADAR*", "Selective – question*",
               "Inclusive (random)", "Inclusive (self-select)")),
  bha_params        = bha_params,
  bha_other_sens    = bha_other_sens,
  nnt_main           = nnt_table_for(main_keys_by_intensity),
  nnt_sens_selective = nnt_table_sens_selective,
  nnt_sens_inclusive = nnt_table_sens_inclusive,
  nnt_pcp            = drop_early(nnt_table_for(pcp_keys_ordered)),
  results_main       = results_table_for(main_keys_by_intensity),
  results_sens_selective = results_table_sens_selective,
  results_sens_inclusive = results_table_sens_inclusive,
  results_pcp            = results_table_for(pcp_keys_ordered),
  timing_burden      = timing_burden_table,
  coverage           = coverage_summary
), file.path(report_dir, "report_data.rds"))

cat("\nWrote", file.path(report_dir, "report_data.rds"), "\n")
