# docs/_s2_figures.R
#
# Draws Figures S2.1-S2.4 of the calibration and validation supplement as PNGs in
# docs/figures/. The supplement runs this itself when it knits (and only when a figure
# is missing or older than its inputs), so there is normally no need to call it by hand.
# To force a redraw, delete docs/figures/ or run, from the project root:
#
#   Rscript docs/_s2_figures.R
#
# Nothing is re-simulated. Every figure is drawn from the saved N = 100,000 benchmarking
# run that analyses/model_development/code/gram_paper1_modelvalidation.Rmd writes, which is
# also the run the numbers in the supplement's text come from. After a recalibration,
# re-run that notebook and point calib_run below at the new file. The check further down
# refuses to draw from a run whose parameters are not the current calibrated ones.
#
# Underscore-prefixed so knitr and rmarkdown::render_site() treat it as a helper, not a
# document to render.

source("model/setup.R")   # opens with rm(list = ls()), so nothing can be assigned above it
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")
source("model/helpers/source_all.R")

calib_run <- "20261112_150020"
sim_file  <- file.path("analyses/model_development/code",
                       paste0("sim_calib_", calib_run, ".RDS"))
if (!file.exists(sim_file)) {
  stop(sim_file, " not found. It is not tracked by git: knit ",
       "analyses/model_development/code/gram_paper1_modelvalidation.Rmd to create it, ",
       "then set calib_run in docs/_s2_figures.R to the run id it prints.")
}
sim_calib <- readRDS(sim_file)

# Refuse a run made under different calibrated parameters: the figures would then
# describe a model the text no longer does.
calib_names <- c("param1", "param2a", "param2b")
ran_with    <- unlist(sim_calib$sim$inputs[calib_names])
current     <- unlist(l.inputs[calib_names])
if (!isTRUE(all.equal(ran_with, current))) {
  stop("The saved run ", calib_run, " was made with ",
       paste(calib_names, "=", ran_with, collapse = ", "),
       ", but the current calibration is ",
       paste(calib_names, "=", current, collapse = ", "),
       ". Re-run gram_paper1_modelvalidation.Rmd and update calib_run.")
}

n_calib <- dim(sim_calib$sim$output)[3]   # only reaches the subtitles, which are dropped

fig_dir <- "docs/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Sized for the 6.5 in text block of a portrait Word page. A ggplot's type scales by
# (display width / render width), so rendering near 6.5 in keeps the type size as set.
save_fig <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), plot = p,
         width = width, height = height, dpi = 300, bg = "white")
}


## S2.1: modeled vs published prevalence of MCI and dementia ####
fig <- compare_prevalence(sim_calib$sim, description = "", n = n_calib)$fig_prev_by_age +
  labs(title = NULL, subtitle = NULL) +
  xlim(NA, 95) + ylim(NA, 0.5)
save_fig(fig, "figS2_1_prevalence_benchmarks", 7.5, 4.6)


## S2.3: predicted cumulative mortality ####
# The helper labels the y axis "Mortality Rate" and adds a caption; what it plots is
# cumulative deaths per 1,000 of the starting cohort, so both are replaced here.
fig <- compare_mortality(sim_calib$sim, description = "", n = n_calib)$cumulative_mort[[1]] +
  labs(title = NULL, subtitle = NULL, caption = NULL,
       y = "Cumulative deaths per 1,000")
save_fig(fig, "figS2_3_mortality_benchmarks", 7.5, 4.6)


## S2.4: time in MCI and dementia by age of onset ####
# Redrawn from compare_reside_time()'s returned data rather than its figure, because the
# reporting version differs from the calibration one in three ways:
#   - the unadjusted "Overall" bar is dropped, leaving only the age-weighted one, which
#     then takes the name "Overall"
#   - empty factor levels sit either side of it, so the two summary bars read as summaries
#     rather than as two more onset bands
#   - no title or subtitle; the caption carries that
# The helper itself is shared with calibration and is left alone.
reside <- compare_reside_time(sim_calib$sim, description = "", n = n_calib)$dat

onset_bins <- setdiff(levels(factor(reside$age_group)),
                      c("Overall", "Overall_adj", "Benchmark"))
x_levels   <- c(onset_bins, " ", "Overall", "  ", "Benchmark")

reside <- reside |>
  subset(age_group != "Overall") |>
  transform(age_group = factor(ifelse(as.character(age_group) == "Overall_adj",
                                      "Overall", as.character(age_group)),
                               levels = x_levels),
            condition = factor(as.character(condition),
                               levels = c("any_dem", "mci", "benchmark_dem", "benchmark_mci")))

fig <- ggplot(reside, aes(x = age_group, y = duration, fill = condition)) +
  geom_col(width = 0.6) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_manual(
    labels = c(any_dem = "Model: Dementia", mci = "Model: MCI",
               benchmark_dem = "Benchmark: Dementia", benchmark_mci = "Benchmark: MCI"),
    values = c(any_dem = "darkgreen", mci = "hotpink",
               benchmark_dem = "lightgreen", benchmark_mci = "pink")) +
  labs(x = "Age of Onset", y = "Duration (years)", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5))
save_fig(fig, "figS2_4_reside_times", 7.5, 4.2)


## S2.2: dementia prevalence stratified by personal characteristics ####
# SEX is stored 1/2; stratify_prevalence_by() expects the 0/1 coding its strat_labels
# are given in. Recoded in place, as in the paper 1 notebook, so this must come after
# the figures above.
sim_calib$sim$output[, "SEX", ] <- sim_calib$sim$output[, "SEX", ] - 1

strata <- list(
  list(var = "SEX",     labels = c("Male", "Female"),                        cutoffs = NULL,
       title = "(A) Sex"),
  list(var = "RACEETH", labels = c("Non-Hispanic White", "Non-Hispanic Black",
                                   "Hispanic", "Other Race"),                cutoffs = NULL,
       title = "(B) Race/Ethnicity"),
  list(var = "INCOME",  labels = c("Low", "Medium", "High"),                 cutoffs = NULL,
       title = "(C) Income Group"),
  list(var = "EDU",     labels = c("Less than HS", "HS + Some College",
                                   "College or More"),        cutoffs = c(0, 12, 16, Inf),
       title = "(D) Educational Attainment"),
  list(var = "MEDBUR",  labels = c("None", "Single Comorbidity",
                                   "Multiple Comorbidity"),    cutoffs = c(0, 1, 2, Inf),
       title = "(E) Medical Burden"),
  list(var = "APOE4",   labels = c("Not Carrier", "Carrier"),                cutoffs = NULL,
       title = "(F) APOE4 Carrier Status")
)

# Each panel's stratum legend sits under it: inset in the top-left corner, as in the
# paper 1 notebook, it covered the curves in four of the six panels at this size. The
# colour legend (dementia vs MCI) is the same in every panel, so it is dropped and the
# caption carries it.
#
# The limits and the type size are added to each panel here rather than to the
# assembled patchwork with `&`: under ggplot2 4.x that operator no longer dispatches on
# a patchwork, which is what the paper 1 notebook's version of this figure used.
panels <- lapply(strata, function(s) {
  stratify_prevalence_by(sim_calib$sim, s$var, strat_labels = s$labels,
                         strat_cutoffs = s$cutoffs, stratify_by_baseline = TRUE)$plot_prev +
    guides(linetype = guide_legend(title = NULL, ncol = 2,
                                   override.aes = list(linewidth = 0.6)),
           color    = "none") +
    theme(legend.position = "bottom",
          legend.margin   = margin(0, 0, 0, 0),
          legend.key.width = unit(1.6, "lines")) +
    labs(title = s$title) +
    xlim(NA, 90) + ylim(NA, 0.4) +
    theme(text        = element_text(size = 10),
          plot.title  = element_text(size = 11, face = "bold"),
          legend.text = element_text(size = 7.5))
})

fig <- patchwork::wrap_plots(panels, ncol = 2)
save_fig(fig, "figS2_2_stratified_prevalence", 7.5, 9.5)

cat("\nWrote Figures S2.1-S2.4 to", fig_dir, "from run", calib_run, "\n")
