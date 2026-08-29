######################################## GRAM-US PAPER 2: BHA TESTING STRATEGIES ########################################

## Setup ####
source("model/setup.R")
source("model/simulation.R")
source("calibration/benchmarking_helpers.R")
source("model/helpers/source_all.R")
source("analyses/testing_strategies/test_performance_helpers.R")
library(flextable)
sample1 <- readRDS("data/acs_data/acs_age50_RACE-revised.RDS")

# Calibrated parameters are applied by model/setup.R, so only the cohort size for this
# analysis is set here.
l.inputs_calibrated <- l.inputs
l.inputs_calibrated[["n.ind"]] <- 100000

# Analysis plan:
# Strategies analyzed will follow the GRAM-ish one-time testing paper
#   (1) reactive testing (with rr 2.0 at a later cycle if previously opted in)
#   (2) selective with eRADAR/EHR-based (with rr 2.0 at a later cycle if previously flagged)
#   (3) inclusive 50% random opt-in (with rr 2.0 at a later cycle if previously opted in)

# Outcomes reported:
#   (a) overall sensitivity/specificity (i.e., strategy level)
#   (b) overall PPV and NPV
#   all at two time points: at the first cycle in which testing is possible, and at end of program (age 80)

# Sensitivity analyses
#   (1) a question-based selective option
#   (2) inclusive with non-random selection (NB, unlike GRAM-ish, total testing not fixed at 50% here)

# PCP follow-up set
#   inclusive and selective, each with an imperfect and a perfect PCP confirming BHA results
#   no PCP follow-up scenario for reactive: PCP already heavily involved in referring to BHA
#                                           so would have no reason to object to a test result


## Scenario registry ####
# key    = short name used in every file name and in the `scenario` column
# config = config file in bha_scenarios/
# label  = facet/table label
# group  = which figure the scenario belongs to
scenario_registry <- tibble::tribble(
  ~key,                   ~config,                                  ~label,                      ~group,
  "u3bhapos_rand50",      "scenario_u3bhapos_rand50_config.R",      "Inclusive",                 "main",
  "s1bhapos_emr",         "scenario_s1bhapos_emr_config.R",         "Selective",                 "main",
  "r1bhapos",             "scenario_r1bhapos_config.R",             "Reactive",                  "main",

  "s1bhapos_question",    "scenario_s1bhapos_question_config.R",    "Selective, question-based", "sens",
  "u3bhapos_nonrand",     "scenario_u3bhapos_nonrand_config.R",     "Inclusive, non-random",     "sens",

  "u3pcppos_rand50",      "scenario_u3pcppos_rand50_config.R",      "Inclusive, imperfect PCP",  "pcp",
  "s1pcppos_emr",         "scenario_s1pcppos_emr_config.R",         "Selective, imperfect PCP",  "pcp",
  "u3pcp_perfect_rand50", "scenario_u3pcp_perfect_rand50_config.R", "Inclusive, perfect PCP",    "pcp",
  "s1pcp_perfect_emr",    "scenario_s1pcp_perfect_emr_config.R",    "Selective, perfect PCP",    "pcp"
)

# Subset this to re-run only part of the set
scenarios_to_run <- scenario_registry$key[1:3]

reg_row    <- function(key) scenario_registry[match(key, scenario_registry$key), ]
labels_for <- function(keys) setNames(reg_row(keys)$label, keys)
keys_in    <- function(...) {
  scenario_registry$key[scenario_registry$group %in% c(...) &
                          scenario_registry$key %in% scenarios_to_run]
}


## Analysis groups and shared constants ####
# Declared here, not inside the figure block, so tables do not depend on plots having run
main_keys           <- keys_in("main")
# Literal rather than keys_in("sens"): each pairs a main-analysis scenario with its
# sensitivity variant, which the registry groups cannot express. Not filtered by
# scenarios_to_run -- the figure loop checks availability instead.
sens_selective_keys <- c("s1bhapos_emr", "s1bhapos_question")
sens_inclusive_keys <- c("u3bhapos_rand50", "u3bhapos_nonrand")
pcp_keys            <- keys_in("pcp")

# A comparison is either complete or it is not produced. Two ways a subset run used
# to go wrong: a group with no scenarios in the run yielded an empty figure, which
# ggsave turned into a blank jpeg and an error that killed the rest of the loop, or a
# 0-row table; and a group missing only SOME of its keys silently rendered a partial
# comparison that looked complete. NB an empty `keys` has nothing "missing" --
# keys_in() returns character(0) when no scenario of that group ran -- so the empty
# case has to be tested separately.
skip_incomplete <- function(keys, what) {
  missing <- setdiff(keys, scenarios_to_run)
  if (!length(keys) || length(missing)) {
    message("Skipping ", what, ": ",
            if (!length(keys)) "no scenarios from this group are in the current run."
            else paste0("not in this run -- ", paste(missing, collapse = ", "),
                        ". Add them to scenarios_to_run, or set run_id to a run that has them."))
    return(TRUE)
  }
  FALSE
}

# flextable() on a skipped (NULL) table errors, so the tail of the script goes
# through this instead.
as_flex <- function(x) if (is.null(x)) invisible(NULL) else flextable(x)

testing_window <- 65:80   # matches age_first_test / age_stop_test in every config
plot_ages      <- 65:80
plot_y_max     <- 75000
end_age        <- 80      # end of follow up period for reporting


## Paths and run id ####
config_dir <- "analyses/testing_strategies/bha_scenarios"
output_dir <- "analyses/testing_strategies/sim_results"
perf_dir   <- "analyses/testing_strategies/test_perf_results"
plot_dir   <- "analyses/testing_strategies/plots"

# To pick up an earlier run, set run_id to that timestamp and skip the "Run scenarios" chunk.
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

sim_file  <- function(key, id = run_id) file.path(output_dir, paste0("scenario_", key, "_sim_", id, ".rds"))
perf_file <- function(key, id = run_id) file.path(perf_dir, paste0(key, "_", id, ".rds"))
stats_file <- function(id = run_id) file.path(perf_dir, paste0("strategy_stats_", id, ".rds"))
plot_file <- function(name, id = run_id) file.path(plot_dir, paste0(name, "_", id, ".jpeg"))


## Run scenarios ####
for (scen in scenarios_to_run) {
  local({
    config_file <- file.path(config_dir, reg_row(scen)$config)
    cat("Running scenario:", scen, "\n")
    config <- load_scenario(config_file, l.inputs_calibrated)
    result <- f.wrap_run(config, microdata = sample1)
    saveRDS(result, file = sim_file(scen))

    rm(config, result)
    invisible(gc())
  })
}


## Test performance data ####
# One pass per scenario produces every derived product that needs the raw array:
# the per-cycle counts the figures and tables run on, the programme size quoted
# in the figure strips, and the yield counts behind the NNT table.
#
#   n_eligible    individuals eligible at ANY point in the window: alive, with a
#                 provider, and no diagnosis as of the previous cycle. The
#                 denominator for programme reach.
#   n_people      individuals tested at least once in the window. The denominator
#                 of the composition row of the counts-and-shares figure. Scoped
#                 to the window rather than the whole run, which is the same thing
#                 only while every config's age_first_test / age_stop_test match
#                 testing_window -- they do today, and no test fires outside it.
#   n_tests       test events in the window. Much the larger number wherever a
#                 strategy re-tests annually.
#   n_identified  individuals who at some point in the window hold a positive
#                 verdict WHILE impaired -- i.e. the strategy got them right.
#                 Distinct people, not the end-of-window TP stock, which loses
#                 people to death and would flatter NNT the later you read it.
#   n_at_mci /    the same people, split on SEV at the cycle they were FIRST
#   n_at_dem      correctly identified: caught at MCI vs caught at dementia.
#   n_early_catch of those, the ones already flagged positive BEFORE they were
#                 impaired -- flagged while healthy OR while in TCI. See the NNT
#                 table for why this one needs care.
#   n_early_at_tci the subset of those flagged during TCI rather than while
#                 outright healthy. Not in any table; carried so the split behind
#                 the NNT footnote can be read off without re-deriving it.
#
# Cached, but only reused if BOTH caches carry every column the reporting code
# reads -- otherwise adding a statistic here, or a column to
# post_processing_outputs(), would silently serve stale results.
stats_cols <- c("scenario", "n_eligible", "n_people", "n_tests", "n_identified",
                "n_at_mci", "n_at_dem", "n_early_catch", "n_early_at_tci")
perf_cols  <- c("tp_mci", "tp_dem", "converted_tp_mci", "converted_tp_dem",
                "fn_mci", "fn_dem", "notest_fn_mci", "notest_fn_dem",
                "not_eligible", "alive")

stats_cached <- file.exists(stats_file()) &&
  all(stats_cols %in% names(readRDS(stats_file()))) &&
  all(file.exists(vapply(scenarios_to_run, perf_file, character(1)))) &&
  all(perf_cols %in% names(readRDS(perf_file(scenarios_to_run[1]))))

strategy_stats <- if (stats_cached) {
  readRDS(stats_file())
} else {
  cyc <- testing_window - 50 + 1

  out <- do.call(rbind, lapply(scenarios_to_run, function(scen) {
    output <- readRDS(sim_file(scen))$output
    saveRDS(post_processing_outputs(output), file = perf_file(scen))

    bha <- output[, "BHA", ]; pcp <- output[, "PCP", ]
    syn <- output[, "SYN", ]; sev <- output[, "SEV", ]
    res <- ifelse(!is.na(pcp) & pcp >= 0, pcp, bha)
    dxr <- apply(res, 2, cummax)          # standing verdict, carried forward

    # Eligible at any point in the window. DX is lagged to match f.update_BHA:
    # modules run BHA before DX, so the testing decision sees last cycle's
    # diagnosis. HCARE is current, which also matches the gate.
    elig <- (output[cyc,     "ALIVE", ] == 1) &
            (output[cyc - 1, "DX", ]    == 0) &
            (output[cyc,     "HCARE", ] == 1)

    bha_w <- bha[cyc, ]
    tp_w  <- ((dxr == 1) & (syn == 1))[cyc, ]   # right about an impaired person
    pos_w <- (dxr == 1)[cyc, ]                  # positive, whatever the truth

    first_tp  <- apply(tp_w,  2, function(x) which(x)[1])
    first_pos <- apply(pos_w, 2, function(x) which(x)[1])
    got    <- which(!is.na(first_tp))
    sev_at <- sev[cyc, ][cbind(first_tp[got], got)]

    # Early catch: the standing positive predates the first cycle the person is
    # both positive and impaired. cummax makes first_pos <= first_tp always, so
    # the strict inequality is exactly "not impaired at the moment of flagging" --
    # which is what makes it cover BOTH ways of being early, flagged while healthy
    # and flagged during TCI. Verified against the arrays: every early catch has
    # SYN 0 or 0.5 at first_pos, every non-early one has SYN 1, and TCI is the
    # larger half of the early group in all three main strategies. The figure bands
    # in post_processing_outputs() and anchored_performance() use the same rule.
    early      <- got[first_pos[got] < first_tp[got]]
    syn_at_pos <- syn[cyc, ][cbind(first_pos[early], early)]

    res_row <- data.frame(
      scenario      = scen,
      n_eligible    = sum(apply(elig, 2, any, na.rm = TRUE)),
      n_people      = sum(apply(bha_w >= 0, 2, any, na.rm = TRUE)),
      n_tests       = sum(bha_w >= 0, na.rm = TRUE),
      n_identified  = length(got),
      n_at_mci      = sum(sev_at == 0),
      n_at_dem      = sum(sev_at >= 1),
      n_early_catch  = length(early),
      n_early_at_tci = sum(syn_at_pos == 0.5))

    rm(output, bha, pcp, syn, sev, res, dxr, elig)
    invisible(gc())
    res_row
  }))

  # Ratios are derived, so they always agree with the counts printed beside them.
  out$per_person   <- out$n_tests / out$n_people
  out$nnt10_any    <- 10 * out$n_tests / out$n_identified
  out$nnt10_at_mci <- 10 * out$n_tests / out$n_at_mci
  out$nnt10_at_dem <- 10 * out$n_tests / out$n_at_dem
  saveRDS(out, stats_file())
  out
}

all_test_data <- do.call(rbind, lapply(scenarios_to_run, function(scen) {
  readRDS(perf_file(scen)) %>% mutate(scenario = scen)
}))

test_data_for <- function(keys) all_test_data %>% filter(scenario %in% keys)


## Reporting helpers ####

# The h / mci / dem split, defined once
status_mix <- function(syn, sev) {
  data.frame(status = case_when(syn < 1  ~ "h",
                                sev == 0 ~ "mci",
                                sev >= 1 ~ "dem")) %>%
    count(status, name = "n") %>%
    mutate(prev = n / sum(n))
}

ci_share <- function(mix) sum(mix$prev[mix$status %in% c("mci", "dem")])

# Composition of the undiagnosed (no DX in previous cycle end), living pool at a given age.
undx_pool <- function(scenario_array, age) {
  cycle <- age - 50 + 1
  d <- data.frame(ALIVE = scenario_array[cycle,     "ALIVE", ],
                  DX    = scenario_array[cycle - 1, "DX", ],
                  SYN   = scenario_array[cycle,     "SYN", ],
                  SEV   = scenario_array[cycle,     "SEV", ]) %>%
    filter(ALIVE == 1, DX == 0)
  status_mix(d$SYN, d$SEV)
}

# Performance at a per-person anchor rather than a common calendar age, so a staggered
# rollout needs no special case. Two anchors are available:
#
#   anchor = "visit"  each person's FIRST CYCLE due: the first cycle they are due
#                     for a test, whether or not one happens. Preferred.
#   anchor = "test"   each person's FIRST ACTUAL TEST, whenever that falls. Kept for
#                     the "who did each strategy actually test?" question.
#
# BHA codes: -9 not due (no provider, prior diagnosis, outside the age window, not this
# person's turn under cohort_split, or already stopped after a positive); -8 due, but
# the concern/uptake gate did not fire; 0/1 tested negative/positive.
#
# The two anchors answer different questions. Under "visit", people who are due but not
# tested stay in the denominator as misses, so the metrics are PROGRAM-level and
# comparable with the year-10 row. Under "test" nobody is untested by construction, so
# notest_tn/notest_fn are 0 and the metrics are TEST-level -- do not put them in the
# same table as the programme-end figures.
#
# Positives keep the paper's early_pos / fp distinction: a positive in someone not yet
# impaired -- healthy or in TCI -- who is impaired later is an early catch, not a plain
# false positive. converted_tp is
# structurally empty -- "positive while healthy, impaired now" cannot happen within one
# cycle -- so it is fixed at 0 and the formulas are predictive_value's with it dropped.
#
# Sensitivity is also reported split by the severity of the impairment being missed:
# SEV 0 is MCI, SEV >= 1 is dementia of any stage. TCI never enters either, because
# impairment here is SYN == 1 and TCI is SYN == 0.5 -- so TCI sits with the healthy,
# which is how every other definition in this analysis treats it. Specificity takes no
# such split: its denominator is the unimpaired, who have no severity.
anchored_performance <- function(scenario_array, anchor = c("visit", "test")) {
  anchor <- match.arg(anchor)
  bha <- scenario_array[, "BHA", ]
  pcp <- scenario_array[, "PCP", ]
  syn <- scenario_array[, "SYN", ]
  sev <- scenario_array[, "SEV", ]
  result <- ifelse(!is.na(pcp) & pcp >= 0, pcp, bha)  # PCP verdict lands in the same cycle

  will_be_impaired <- colSums(syn == 1, na.rm = TRUE) > 0   # ever impaired, whole horizon

  # -8 counts as an anchor for "visit", but not for "test"
  cutoff <- if (anchor == "visit") -8 else 0
  first_visit <- apply(bha, 2, function(x) which(x >= cutoff)[1])
  visited <- which(!is.na(first_visit))
  idx <- cbind(first_visit[visited], visited)

  res      <- result[idx]
  impaired <- syn[idx] == 1
  early    <- will_be_impaired[visited]
  tested   <- res >= 0
  at_mci   <- impaired & sev[idx] == 0    # SEV only means anything where impaired,
  at_dem   <- impaired & sev[idx] >= 1    # so both are gated on it

  perf <- data.frame(
    n_visited = length(visited),
    n_tested  = sum(tested),
    tp = sum(res ==  1 &  impaired),
    fp = sum(res ==  1 & !impaired & !early),   # positive, never impaired
    early_pos = sum(res == 1 & !impaired & early),   # positive, impaired later
    converted_tp = 0,                                # impossible within one cycle
    tn = sum(res ==  0 & !impaired),
    fn = sum(res ==  0 &  impaired),
    notest_tn = sum(res == -8 & !impaired),
    notest_fn = sum(res == -8 &  impaired),
    tp_mci        = sum(res ==  1 & at_mci),
    fn_mci        = sum(res ==  0 & at_mci),
    notest_fn_mci = sum(res == -8 & at_mci),
    tp_dem        = sum(res ==  1 & at_dem),
    fn_dem        = sum(res ==  0 & at_dem),
    notest_fn_dem = sum(res == -8 & at_dem)
  ) %>%
    mutate(sens = tp / (tp + fn + notest_fn),
           sens_mci = tp_mci / (tp_mci + fn_mci + notest_fn_mci),
           sens_dem = tp_dem / (tp_dem + fn_dem + notest_fn_dem),
           spec = (tn + notest_tn) / (tn + notest_tn + early_pos + fp),
           ppv  = tp / (tp + early_pos + fp),
           npv  = (tn + notest_tn) / (tn + notest_tn + fn + notest_fn),
           acc  = (tp + tn + notest_tn) / (tp + early_pos + fp + tn + fn + notest_tn + notest_fn))

  list(
    anchor     = anchor,
    perf       = perf,
    tested_mix = status_mix(syn[idx][tested], sev[idx][tested]), # cognitive mix among those who actually receive the test
    visit_age  = table(scenario_array[, "AGE", ][idx])
  )
}

first_visit_performance <- function(scenario_array) anchored_performance(scenario_array, "visit")
first_test_performance  <- function(scenario_array) anchored_performance(scenario_array, "test")


## Figures ####
theme_paper2 <- theme(
  text            = element_text(size = 18),
  axis.title      = element_text(size = 20, face = "bold"),
  axis.text       = element_text(size = 16),
  legend.text     = element_text(size = 16),
  legend.title    = element_text(size = 18, face = "bold"),
  strip.text      = element_text(size = 20, face = "bold"),
  plot.title      = element_text(size = 24, face = "bold"),
  legend.key.size = unit(1.2, "cm")
)

# main_keys follows registry order (Inclusive first); the counts+shares figure
# reads least to most intensive, so it takes its own explicit ordering.
main_keys_by_intensity <- c("r1bhapos", "s1bhapos_emr", "u3bhapos_rand50")

figure_specs <- list(
  list(name = "counts-and-shares",           keys = main_keys_by_intensity,
                                                                 type = "combined", width = 17, height = 12.5),
  list(name = "no-early-positives",          keys = main_keys,           type = "results", width = 14),
  list(name = "testers",                     keys = main_keys,           type = "testers", width = 14),
  list(name = "sens-selective",              keys = sens_selective_keys, type = "results", width = 11),
  list(name = "sens-inclusive",              keys = sens_inclusive_keys, type = "results", width = 11),
  list(name = "WITH-PCP-no-early-positives", keys = pcp_keys,            type = "results", width = 14),
  list(name = "WITH-PCP-testers",            keys = pcp_keys,            type = "testers", width = 14)
)

# Draws and saves every panel; returns them so any one can be viewed, e.g.
#   figures[["sens-selective"]]
figures <- setNames(lapply(figure_specs, function(spec) {
  if (skip_incomplete(spec$keys, paste0("figure '", spec$name, "'"))) return(NULL)

  d    <- test_data_for(spec$keys)
  labs <- labels_for(spec$keys)

  p <- switch(spec$type,
    results  = plot_test_results(d, ages = plot_ages, show_early_pos = FALSE,
                                 scenario_names = labs, y_max = plot_y_max),
    testers  = plot_testers(d, ages = plot_ages, scenario_names = labs),
    combined = plot_counts_and_shares(d, scenario_names = labs, ages = plot_ages,
                                      strategy_stats = strategy_stats[strategy_stats$scenario %in% spec$keys, ],
                                      base_size = 20, label_size = 4.6))

  # plot_counts_and_shares carries its own theme, sized via base_size. Adding
  # theme_paper2 on top would override the bottom row's small grey stats strip.
  themed <- if (spec$type == "combined") p else p + theme_paper2

  ggsave(plot_file(spec$name), plot = themed,
         height = spec$height %||% 10, width = spec$width, dpi = 300)
  p
}), vapply(figure_specs, `[[`, character(1), "name"))

figures[["counts-and-shares"]]
figures[["no-early-positives"]]
names(figures)
figures[["testers"]] + ylim(NA, 75000)

## Reporting: methods ####
# Table 1: testing likelihood by strategy and cognitive state.
# Read from the configs the run actually used, so the table cannot drift from the
# scenarios. This replaces a hand-assembled version that quoted the question-based
# probabilities on the selective row while the main analysis ran the eRADAR ones.
testing_likelihood_for <- function(keys) {
  do.call(rbind, lapply(keys, function(k) {
    scen <- load_scenario(file.path(config_dir, reg_row(k)$config), l.inputs_calibrated)[["scenario"]]
    m <- scen[["probs_select"]][1, ]
    data.frame(Strategy = reg_row(k)$label,
               Healthy = m$h, MCI = m$mci, Dementia = m$dem,
               `RR if concern last cycle` = scen[["rr.select_prior"]],
               check.names = FALSE)
  }))
}

tab1 <- flextable(testing_likelihood_for(c(main_keys, sens_selective_keys[-1],
                                           sens_inclusive_keys[-1])))
tab1

# Test characteristics quoted in the methods. BHA_GS comes from model/test_properties.R;
# rr.select_prior is per-scenario now and appears in tab1 above.
bha_test_params <- list(
  sensitivity = BHA_GS$sens,
  specificity = BHA_GS$spec
)
bha_test_params


## Reporting: results ####
# Composition of the undiagnosed pool at the ages quoted in the text. Read from one
# reference scenario so the quoted figures describe a single cohort.
denom_ref_key <- "u3bhapos_rand50"
denom_ref     <- readRDS(sim_file(denom_ref_key))$output
prev_in_undx  <- setNames(lapply(c(65, 67), function(a) undx_pool(denom_ref, a)),
                          as.character(c(65, 67)))
rm(denom_ref); invisible(gc())

prev_in_undx

# One pass over each scenario's array: everything that needs the raw output
scenario_reports <- setNames(lapply(scenarios_to_run, function(scen) {
  arr <- readRDS(sim_file(scen))$output
  out <- list(
    first_visit = first_visit_performance(arr),   # reported
    first_test  = first_test_performance(arr)     # alternative, if asked for
  )
  rm(arr); invisible(gc())
  out
}), scenarios_to_run)

# Coverage over the life of the program, and who the first round reached
coverage_summary <- do.call(rbind, lapply(scenarios_to_run, function(scen) {
  r  <- scenario_reports[[scen]]
  fv <- r$first_visit
  st <- strategy_stats[match(scen, strategy_stats$scenario), ]
  data.frame(
    scenario        = scen,
    strategy        = reg_row(scen)$label,
    n_ever_eligible = st$n_eligible,
    # Both terms span the whole testing window: of everyone the programme could have
    # reached, what share ever received a test
    n_ever_tested       = st$n_people,
    pct_eligible_tested = st$n_people / st$n_eligible,
    # The first round only
    n_visited           = fv$perf$n_visited,
    n_tested_first      = fv$perf$n_tested,
    ci_prev_first       = ci_share(fv$tested_mix)
  )
}))

coverage_summary

# Cognitive status of those tested, by strategy. anchor = "first_visit" is the reported
# version (the first round only); "first_test" is everyone ever tested, read at whenever
# their own first test fell -- kept in case a reviewer asks for it.
tested_mix_for <- function(anchor = c("first_visit", "first_test")) {
  anchor <- match.arg(anchor)
  do.call(rbind, lapply(scenarios_to_run, function(scen) {
    scenario_reports[[scen]][[anchor]]$tested_mix %>%
      mutate(strategy = reg_row(scen)$label, .before = 1)
  }))
}

first_visit_mix <- tested_mix_for("first_visit")
first_visit_mix


# Counts behind the figures. NB this folds early_pos into fp and converted_tp into tp,
# so the early-catch split is deliberately collapsed here; use scenario_reports for it.
counts_table <- function(keys, ages = c(65, 67, 75)) {
  test_data_for(keys) %>%
    filter(age %in% ages) %>%
    mutate(age = age,
           scenario = scenario,
           tp = tp + converted_tp,
           fp = fp + early_pos,
           tn = tn + notest_tn,
           fn = fn + notest_fn,
           dead = death,
           .keep = "none") %>%
    flextable()
}

counts_table(main_keys)


# Program performance by calendar age, used for the year-10 row
predictive_value <- all_test_data %>%
  mutate(ppv = (tp + converted_tp) / (tp + early_pos + converted_tp + fp),
         npv = ((tn + notest_tn) / (tn + fn + notest_tn + notest_fn)),
         sens = (tp + converted_tp) / (tp + converted_tp + fn + notest_fn),
         sens_mci = (tp_mci + converted_tp_mci) /
                      (tp_mci + converted_tp_mci + fn_mci + notest_fn_mci),
         sens_dem = (tp_dem + converted_tp_dem) /
                      (tp_dem + converted_tp_dem + fn_dem + notest_fn_dem),
         spec = (tn + notest_tn) / (tn + notest_tn + early_pos + fp),
         acc = (tp + converted_tp + tn + notest_tn) / (tp + converted_tp + tn + notest_tn + fp + early_pos + fn + notest_fn))

# `eligible_tested` is supplied by the caller rather than read off `d`, because its
# denominator is row-specific -- see results_table_for.
as_pct_row <- function(d, time, strategy, eligible_tested) {
  pct <- function(x) round(x * 100, digits = 1)
  data.frame(Time                    = time,
             Strategy                = strategy,
             `Eligible tested`       = pct(eligible_tested),
             Sensitivity             = pct(d$sens),
             `Sensitivity, MCI`      = pct(d$sens_mci),
             `Sensitivity, dementia` = pct(d$sens_dem),
             Specificity             = pct(d$spec),
             PPV                     = pct(d$ppv),
             NPV                     = pct(d$npv),
             Accuracy                = pct(d$acc),
             check.names = FALSE)
}

# "First Cycle" is read at each person's own first cycle due for a test, so no calendar
# age is involved. "By Program End" is the calendar snapshot at end_age. Both are
# program-level: eligible people who were never tested count as misses in each.
# anchor = "first_test" swaps in the test-level first row; see anchored_performance for
# why that row is NOT comparable with the program-end one.
#
# "Eligible tested" is the share of eligible people the strategy actually tested. Its
# denominator is MATCHED TO THE ROW, so the column reads consistently as reach-so-far:
# on the First Cycle row, of those due a test that cycle, the share tested; on the
# Program End row, of everyone ever eligible across the whole window, the share ever
# tested. Under anchor = "first_test" the first row is 100% by construction -- nobody
# is untested there -- which is one more reason that row does not belong beside the
# program-end one.
results_table_for <- function(keys, at_end_followup = end_age,
                              anchor = c("first_visit", "first_test")) {
  if (skip_incomplete(keys, "results table")) return(NULL)
  anchor <- match.arg(anchor)
  first_label <- if (anchor == "first_visit") "First Cycle" else "First Test"

  first <- do.call(rbind, lapply(keys, function(scen) {
    fv <- scenario_reports[[scen]][[anchor]]$perf
    as_pct_row(fv, first_label, reg_row(scen)$label,
               eligible_tested = fv$n_tested / fv$n_visited)
  }))

  later <- do.call(rbind, lapply(keys, function(scen) {
    st <- strategy_stats[match(scen, strategy_stats$scenario), ]
    as_pct_row(predictive_value %>% filter(scenario == scen, age == at_end_followup),
               "By Program End", reg_row(scen)$label,
               eligible_tested = st$n_people / st$n_eligible)
  }))

  rbind(first, later)
}

# Number needed to test: how many BHAs a strategy runs per 10 people it correctly
# identifies. All three NNT columns keep the SAME numerator -- every test the
# strategy runs -- because no test can be aimed at MCI alone; only the denominator
# changes. So the MCI and dementia columns read as "tests needed to find 10 people
# at that stage", and they do not average to the overall column.
#
# "Correctly identified" is the first cycle a person holds a positive verdict
# while impaired, and the MCI/dementia split is their SEV at that cycle -- caught
# at MCI, not has MCI now.
#
# The last column counts EARLY catches, and it counts both ways of being early:
# flagged while outright healthy and later converting, and flagged during TCI and
# later converting. Both fall out of the same test -- the standing positive
# predates the first impaired-and-positive cycle -- because impairment here is
# SYN == 1 and TCI is SYN == 0.5, so a TCI-cycle positive is not yet a true
# positive. Checked against the arrays: every early catch has SYN 0 or 0.5 when
# flagged, every non-early one has SYN 1, and TCI is roughly two thirds of the
# early group in each main strategy (strategy_stats$n_early_at_tci has the split).
#
# UNRESOLVED, for the team: anyone flagged before they were impaired necessarily
# enters the identified state at the moment they convert, which is by definition
# MCI -- so an early catch can ONLY ever land in the MCI column, never in
# dementia. Whether those belong in "caught at MCI", in a column of their own, or
# outside the NNT denominator altogether is a judgement about what the paper is
# claiming, not a coding question. Shown as a separate count for now so the
# choice is visible rather than buried.
nnt_table_for <- function(keys) {
  if (skip_incomplete(keys, "NNT table")) return(NULL)
  d <- strategy_stats[match(keys, strategy_stats$scenario), ]
  data.frame(
    Strategy                    = reg_row(keys)$label,
    Tests                       = d$n_tests,
    `Correctly identified`      = d$n_identified,
    `Tests per 10 identified`   = round(d$nnt10_any, 1),
    `Caught at MCI`             = d$n_at_mci,
    `Tests per 10 at MCI`       = round(d$nnt10_at_mci, 1),
    `Caught at dementia`        = d$n_at_dem,
    `Tests per 10 at dementia`  = round(d$nnt10_at_dem, 1),
    `of MCI: flagged before impairment` = d$n_early_catch,
    check.names = FALSE)
}

nnt_table <- nnt_table_for(main_keys)
as_flex(nnt_table)

# Same table for the other groups, if the call runs long
nnt_table_sens_selective <- nnt_table_for(sens_selective_keys)
nnt_table_sens_inclusive <- nnt_table_for(sens_inclusive_keys)
nnt_table_pcp            <- nnt_table_for(pcp_keys)


results_table <- results_table_for(main_keys)
as_flex(results_table)

results_table_sens_selective <- results_table_for(sens_selective_keys)
as_flex(results_table_sens_selective)

results_table_sens_inclusive <- results_table_for(sens_inclusive_keys)
as_flex(results_table_sens_inclusive)

results_table_pcp <- results_table_for(pcp_keys)
as_flex(results_table_pcp)

