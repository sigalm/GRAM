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
#   inclusive and selective, each with a perfect PCP confirming BHA results (reported), and
#   an imperfect one (registry group "pcp_imperfect", not reported)
#   reactive has the perfect PCP only, for completeness: the PCP is already heavily involved
#                                      in referring to BHA, so it is the arm least likely to
#                                      see a test result overturned


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

  "u3pcp_perfect_rand50", "scenario_u3pcp_perfect_rand50_config.R", "Inclusive, perfect PCP",    "pcp",
  "s1pcp_perfect_emr",    "scenario_s1pcp_perfect_emr_config.R",    "Selective, perfect PCP",    "pcp",
  "r1pcp_perfect",        "scenario_r1pcp_perfect_config.R",        "Reactive, perfect PCP",     "pcp",

  # Not reported: the perfect PCP is the bound on how much confirmation can undo, and
  # the imperfect one adds a second layer of assumptions without a new question.
  # Kept runnable by adding the group to scenarios_to_run.
  "u3pcppos_rand50",      "scenario_u3pcppos_rand50_config.R",      "Inclusive, imperfect PCP",  "pcp_imperfect",
  "s1pcppos_emr",         "scenario_s1pcppos_emr_config.R",         "Selective, imperfect PCP",  "pcp_imperfect"
)

# Subset this to re-run only part of the set. options(gram.scenarios = ...) overrides it,
# for the same reason the run options below are options -- see "Paths and run id".
scenarios_to_run <- getOption("gram.scenarios",
                              scenario_registry$key[scenario_registry$group %in% c("main")])

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

# To pick up an earlier run without re-simulating it, set these BEFORE sourcing this
# file and every table and figure is rebuilt from that run's saved outputs:
#
#   options(gram.run_id = "<timestamp>", gram.run_scenarios = FALSE)
#   options(gram.make_figures = FALSE)   # additionally skip redrawing the jpegs
#
# Options rather than plain variables because model/setup.R, sourced at the top of this
# file, opens with rm(list = ls()) -- anything the caller assigns beforehand is gone by
# the time this line runs. Interactively, assigning run_id after the setup block and
# skipping the "Run scenarios" chunk by hand does the same thing. Every default is the
# interactive behaviour: fresh timestamp, scenarios run, figures written.
run_id        <- getOption("gram.run_id", format(Sys.time(), "%Y%m%d_%H%M%S"))
run_scenarios <- getOption("gram.run_scenarios", TRUE)
make_figures  <- getOption("gram.make_figures", TRUE)

# Each run gets its own subfolder in each of the three output folders. The run id stays
# in the file name too, so a file copied out of its folder still says which run made it.
sim_file  <- function(key, id = run_id) file.path(output_dir, id, paste0("scenario_", key, "_sim_", id, ".rds"))
perf_file <- function(key, id = run_id) file.path(perf_dir, id, paste0(key, "_", id, ".rds"))
stats_file <- function(id = run_id) file.path(perf_dir, id, paste0("strategy_stats_", id, ".rds"))
plot_file <- function(name, id = run_id) file.path(plot_dir, id, paste0(name, "_", id, ".jpeg"))

# Wrap a path that is about to be written: creates its run folder on first write, so
# reading an earlier run never leaves an empty folder behind.
for_writing <- function(path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  path
}


## Run scenarios ####
for (scen in if (run_scenarios) scenarios_to_run else character(0)) {
  local({
    config_file <- file.path(config_dir, reg_row(scen)$config)
    cat("Running scenario:", scen, "\n")
    config <- load_scenario(config_file, l.inputs_calibrated)
    result <- f.wrap_run(config, microdata = sample1)
    saveRDS(result, file = for_writing(sim_file(scen)))

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
#   n_at_mci /    of those, the ones already impaired when they tested positive,
#   n_at_dem      split on SEV at that test: caught at MCI vs caught at dementia.
#   n_early_catch the rest: flagged positive BEFORE they were impaired -- while
#                 healthy OR while in TCI. The three groups sum to n_identified.
#   n_early_at_tci the subset of those flagged during TCI rather than while
#                 outright healthy. Not in any table; carried so the split behind
#                 the NNT footnote can be read off without re-deriving it.
#   nnt10_*_lo /  95% bootstrap interval of the NNT table's three ratios, at MCI,
#   nnt10_*_hi    at dementia, and correctly detected (MCI + dementia). See nnt_boot().
#
# Cached, but only reused if BOTH caches carry every column the reporting code
# reads -- otherwise adding a statistic here, or a column to
# post_processing_outputs(), would silently serve stale results.
stats_cols <- c("scenario", "n_eligible", "n_people", "n_tests", "n_identified",
                "n_at_mci", "n_at_dem", "n_early_catch", "n_early_at_tci",
                "nnt10_early", "nnt10_detected",
                "nnt10_at_mci_lo", "nnt10_at_mci_hi", "nnt10_at_dem_lo", "nnt10_at_dem_hi",
                "nnt10_detected_lo", "nnt10_detected_hi")
perf_cols  <- c("tp_mci", "tp_dem", "converted_tp_mci", "converted_tp_dem",
                "fn_mci", "fn_dem", "notest_fn_mci", "notest_fn_dem",
                "notest_fn_dx", "not_eligible", "alive",
                "identified", "identified_mci", "identified_dem", "identified_early",
                "new_tp", "new_tp_converted", "new_pos", "new_pos_impaired")

stats_cached <- file.exists(stats_file()) &&
  all(stats_cols %in% names(readRDS(stats_file()))) &&
  all(file.exists(vapply(scenarios_to_run, perf_file, character(1)))) &&
  all(perf_cols %in% names(readRDS(perf_file(scenarios_to_run[1]))))

# Monte Carlo intervals for the NNT ratios, by bootstrapping people rather than tests:
# one person takes many tests, and their tests and their detection go together. Each
# replicate reweights the cohort by how many times each person is drawn, so it is one
# crossprod rather than a copy of the cohort. The seed is fixed and the same for every
# scenario, and the arms share their simulated people, so every arm is resampled the
# same way. Like the results table's Wilson intervals, these carry Monte Carlo error
# only, not uncertainty in any input. Returns a 2 x k matrix, rows the 2.5th and 97.5th
# percentiles, one column per column of `detected`.
nnt_boot <- function(tests, detected, n_boot = 1000, seed = 20261029) {
  n <- length(tests)
  m <- cbind(tests, detected)
  sums <- withr::with_seed(seed, t(vapply(seq_len(n_boot), function(b)
    drop(crossprod(tabulate(sample.int(n, n, replace = TRUE), n), m)), numeric(ncol(m)))))
  ratios <- 10 * sums[, 1] / sums[, -1, drop = FALSE]
  apply(ratios, 2, quantile, probs = c(0.025, 0.975), names = FALSE)
}

strategy_stats <- if (stats_cached) {
  readRDS(stats_file())
} else {
  cyc <- testing_window - 50 + 1

  out <- do.call(rbind, lapply(scenarios_to_run, function(scen) {
    output <- readRDS(sim_file(scen))$output
    perf <- post_processing_outputs(output, age_first_test = min(testing_window))
    saveRDS(perf, file = for_writing(perf_file(scen)))

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
    is_early   <- first_pos[got] < first_tp[got]
    early      <- got[is_early]
    syn_at_pos <- syn[cyc, ][cbind(first_pos[early], early)]

    # Per person, for the bootstrap: tests taken, and detected at MCI / at dementia
    at_mci <- at_dem <- logical(ncol(bha_w))
    at_mci[got[!is_early & sev_at == 0]] <- TRUE
    at_dem[got[!is_early & sev_at >= 1]] <- TRUE
    ci <- nnt_boot(colSums(bha_w >= 0, na.rm = TRUE),
                   cbind(at_mci, at_dem, at_mci | at_dem))

    # For everyone else first_tp IS the positive test, so sev_at is the stage they
    # were at when they tested positive
    res_row <- data.frame(
      scenario      = scen,
      n_eligible    = sum(apply(elig, 2, any, na.rm = TRUE)),
      n_people      = sum(apply(bha_w >= 0, 2, any, na.rm = TRUE)),
      n_tests       = sum(bha_w >= 0, na.rm = TRUE),
      n_identified  = length(got),
      n_at_mci      = sum(sev_at[!is_early] == 0),
      n_at_dem      = sum(sev_at[!is_early] >= 1),
      n_early_catch  = length(early),
      n_early_at_tci = sum(syn_at_pos == 0.5),
      nnt10_at_mci_lo   = ci[1, 1], nnt10_at_mci_hi   = ci[2, 1],
      nnt10_at_dem_lo   = ci[1, 2], nnt10_at_dem_hi   = ci[2, 2],
      nnt10_detected_lo = ci[1, 3], nnt10_detected_hi = ci[2, 3])

    # The identified-by-stage figure reads the same counts per age from
    # post_processing_outputs(), computed separately; they must meet at window end.
    at_end <- perf[perf$age == max(testing_window), ]
    stopifnot(at_end$identified       == res_row$n_identified,
              at_end$identified_mci   == res_row$n_at_mci,
              at_end$identified_dem   == res_row$n_at_dem,
              at_end$identified_early == res_row$n_early_catch)

    rm(output, perf, bha, pcp, syn, sev, res, dxr, elig, at_mci, at_dem)
    invisible(gc())
    res_row
  }))

  # Ratios are derived, so they always agree with the counts printed beside them.
  out$per_person   <- out$n_tests / out$n_people
  out$nnt10_any    <- 10 * out$n_tests / out$n_identified
  out$nnt10_at_mci <- 10 * out$n_tests / out$n_at_mci
  out$nnt10_at_dem <- 10 * out$n_tests / out$n_at_dem
  out$nnt10_early  <- 10 * out$n_tests / out$n_early_catch
  out$nnt10_detected <- 10 * out$n_tests / (out$n_at_mci + out$n_at_dem)
  saveRDS(out, for_writing(stats_file()))
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
# cycle -- so it is fixed at 0. The metrics themselves are computed from these counts by
# metric_counts(), the same formulas as the program-end row.
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
  )

  list(
    anchor     = anchor,
    perf       = perf,
    tested_mix = status_mix(syn[idx][tested], sev[idx][tested]), # cognitive mix among those who actually receive the test
    visit_age  = table(scenario_array[, "AGE", ][idx])
  )
}

first_visit_performance <- function(scenario_array) anchored_performance(scenario_array, "visit")
first_test_performance  <- function(scenario_array) anchored_performance(scenario_array, "test")


# Time from impairment onset to identification, per person. The cohort is INCIDENT
# cases: people whose impairment (SYN 1) begins inside the testing window while they
# are in the program -- provider, no clinical diagnosis as of the previous cycle.
# People already impaired at 65 are left out, because their onset predates the
# program and the clock would start before it could act.
#
# Each person lands in exactly one outcome by end of window:
#   early       positive before onset (healthy or TCI). A false positive until they
#               convert, a true positive from then on -- as everywhere else -- but
#               kept out of the lag, where it would read as caught on day one.
#   caught      first positive while impaired, at lag = years after onset
#               (0 = the onset year itself)
#   clinic      clinically diagnosed before the program caught them: a miss
#   neither     not caught by end of window, and not diagnosed -- died, or still
#               undetected
# Returns one row per incident case, so the table can summarize it however needed.
onset_timing <- function(scenario_array, window = testing_window) {
  bha <- scenario_array[, "BHA", ]; pcp <- scenario_array[, "PCP", ]
  syn <- scenario_array[, "SYN", ]; age <- scenario_array[, "AGE", ]
  dx  <- scenario_array[, "DX", ];  hcare <- scenario_array[, "HCARE", ]
  res <- ifelse(!is.na(pcp) & pcp >= 0, pcp, bha)   # PCP verdict lands in the same cycle
  pos <- apply(res, 2, cummax) == 1                 # standing positive verdict

  first_of <- function(m) apply(m, 2, function(x) which(x)[1])
  onset     <- first_of(syn == 1)
  first_pos <- first_of(pos)
  first_tp  <- first_of(pos & syn == 1)
  first_dx  <- first_of(dx == 1)
  end_cyc   <- max(window) - 50 + 1

  j <- which(!is.na(onset) & onset > 1)
  j <- j[age[cbind(onset[j], j)] %in% window &
         dx[cbind(onset[j] - 1, j)] == 0 &
         hcare[cbind(onset[j], j)] == 1]

  early  <- !is.na(first_pos[j]) & first_pos[j] < onset[j]
  caught <- !early & !is.na(first_tp[j]) & first_tp[j] <= end_cyc
  clinic <- !early & !caught & !is.na(first_dx[j]) & first_dx[j] <= end_cyc

  data.frame(
    onset_age = age[cbind(onset[j], j)],
    outcome   = ifelse(early, "early", ifelse(caught, "caught", ifelse(clinic, "clinic", "neither"))),
    lag       = ifelse(caught, first_tp[j] - onset[j], NA),
    # Early catches: years they carried a positive before they were impaired -- time
    # spent labeled while unimpaired. Bounded by the window, since nobody is tested
    # before age_first_test.
    lead      = ifelse(early, onset[j] - first_pos[j], NA)
  )
}


# Years each person spent holding a positive result while unimpaired (healthy or in
# TCI), counted over the testing window -- through end of program, though the label
# outlasts it, since the stop rule means no one is retested to clear it. Both kinds
# of false positive: early catches, whose label stops counting when they convert,
# and people never impaired, whose label never does. One row per person ever
# labeled. The yearly stock behind it is fp + early_pos in post_processing_outputs().
labeled_unimpaired <- function(scenario_array, ages = testing_window) {
  bha <- scenario_array[, "BHA", ]; pcp <- scenario_array[, "PCP", ]
  syn <- scenario_array[, "SYN", ]
  res <- ifelse(!is.na(pcp) & pcp >= 0, pcp, bha)
  pos <- apply(res, 2, cummax) == 1

  # SYN is NA after death, so the dead stop counting
  years <- colSums((pos & syn < 1)[ages - 50 + 1, , drop = FALSE], na.rm = TRUE)
  later <- colSums(syn == 1, na.rm = TRUE) > 0   # impaired at any age, as for early_pos
  data.frame(years = years[years > 0], later_impaired = later[years > 0])
}


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
  list(name = "annual-performance",          keys = main_keys_by_intensity,
                                                                 type = "annual",  width = 15, height = 6.5),
  list(name = "identified-by-stage",         keys = main_keys_by_intensity,
                                                                 type = "stages",  width = 17, height = 7),
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
    stages   = plot_identified_by_stage(d, scenario_names = labs, ages = plot_ages,
                                        base_size = 20, label_size = 4.6),
    annual   = plot_annual_performance(d, scenario_names = labs, ages = plot_ages,
                                       base_size = 20, label_size = 5.5),
    combined = plot_counts_and_shares(d, scenario_names = labs, ages = plot_ages,
                                      strategy_stats = strategy_stats[strategy_stats$scenario %in% spec$keys, ],
                                      base_size = 20, label_size = 4.6))

  # plot_counts_and_shares and plot_identified_by_stage carry their own theme,
  # sized via base_size. Adding theme_paper2 on top would override the bottom row's
  # small grey stats strip, and bring back a legend the line-end labels replace.
  themed <- if (spec$type %in% c("combined", "stages", "annual")) p else p + theme_paper2

  if (make_figures)
    ggsave(for_writing(plot_file(spec$name)), plot = themed,
           height = spec$height %||% 10, width = spec$width, dpi = 300)
  p
}), vapply(figure_specs, `[[`, character(1), "name"))

figures[["counts-and-shares"]]
figures[["no-early-positives"]]
figures[["identified-by-stage"]]
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
    # P(tested) at a first offer: selection times uptake. Where a scenario has no
    # probs_accept every offer is taken up, and probs_select is already P(tested).
    m <- scen[["probs_select"]][1, c("h", "mci", "dem")]
    if (!is.null(scen[["probs_accept"]])) m <- m * scen[["probs_accept"]][1, c("h", "mci", "dem")]
    data.frame(Strategy = reg_row(k)$label,
               Healthy = m$h, MCI = m$mci, Dementia = m$dem,
               `Redrawn each year` = if (isTRUE(scen[["select_persists"]])) "No" else "Yes",
               check.names = FALSE)
  }))
}

# Testing uptake inputs, the two parts of testing_likelihood_for() shown apart: offered
# (probs_select), completed if offered (probs_accept, 1 where unset), and completed after
# refusing the last offer (p.accept_after_decline; N/A where no one can refuse). Measures
# down, strategies across; `labels` names the columns, in the order of `keys`. Block is
# for the Word report to group rows on. Read from the configs, like the table above.
uptake_inputs_for <- function(keys, labels) {
  fmt <- function(p) ifelse(p < 0.01, sprintf("%.1f%%", 100 * p), sprintf("%.0f%%", 100 * p))
  cols <- lapply(keys, function(k) {
    scen   <- load_scenario(file.path(config_dir, reg_row(k)$config), l.inputs_calibrated)[["scenario"]]
    sel    <- unlist(scen[["probs_select"]][1, c("h", "mci", "dem")])
    acc    <- if (is.null(scen[["probs_accept"]])) c(h = 1, mci = 1, dem = 1)
              else unlist(scen[["probs_accept"]][1, c("h", "mci", "dem")])
    refuse <- scen[["p.accept_after_decline"]]
    c(fmt(sel), fmt(acc), if (is.null(refuse) || all(acc == 1)) "N/A" else fmt(refuse))
  })
  out <- data.frame(
    Block   = c(rep("Probability of being offered a test", 3),
                rep("Probability of completing a test if offered", 4)),
    Measure = c("Healthy", "MCI", "Dementia", "Healthy", "MCI", "Dementia",
                "Refused last cycle, any cognitive status"),
    check.names = FALSE)
  out[labels] <- cols
  out
}

tab1 <- flextable(testing_likelihood_for(c(main_keys, sens_selective_keys[-1],
                                           sens_inclusive_keys[-1])))
tab1

# Test characteristics quoted in the methods. BHA_GS comes from model/test_properties.R;
# how selection repeats is per-scenario and appears in tab1 above.
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
    first_test  = first_test_performance(arr),    # alternative, if asked for
    timing      = onset_timing(arr),
    labeled     = labeled_unimpaired(arr)
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
counts_table <- function(keys, ages = c(65, 67, 80)) {
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


# Numerator and denominator of every reported metric, from one row of counts: either a
# first-cycle perf (anchored_performance) or one age of all_test_data. The single place
# these formulas live. converted_tp_mci / _dem exist only in all_test_data; the
# first-cycle row has no conversions, so they count as 0 there.
metric_counts <- function(d) {
  with_conv <- function(col) d[[col]] + (d[[paste0("converted_", col)]] %||% 0)
  tp     <- with_conv("tp")
  tp_mci <- with_conv("tp_mci")
  tp_dem <- with_conv("tp_dem")
  neg    <- d$tn + d$notest_tn   # correct negatives, tested or not
  list(sens     = c(x = tp,       n = tp + d$fn + d$notest_fn),
       sens_mci = c(x = tp_mci,   n = tp_mci + d$fn_mci + d$notest_fn_mci),
       sens_dem = c(x = tp_dem,   n = tp_dem + d$fn_dem + d$notest_fn_dem),
       spec     = c(x = neg,      n = neg + d$early_pos + d$fp),
       ppv      = c(x = tp,       n = tp + d$early_pos + d$fp),
       npv      = c(x = neg,      n = neg + d$fn + d$notest_fn),
       acc      = c(x = tp + neg, n = tp + neg + d$early_pos + d$fp + d$fn + d$notest_fn))
}

# Percent with a 95% Wilson interval, as one cell: "50.4 (44.2-56.6)". Wilson rather than
# Wald because the first-cycle dementia cells rest on ~250 people and sit near 1 for PPV
# and specificity, where Wald runs past 100%.
#
# These intervals carry Monte Carlo error only -- would the cell move with a different
# draw of the same cohort size -- not uncertainty in any input parameter. And they are
# not for comparing strategies: the arms share random numbers, so the same simulated
# people sit behind every row, and a difference between arms is far more precise than
# two overlapping intervals suggest.
pct_ci <- function(x, n, z = qnorm(0.975)) {
  p      <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  sprintf("%.1f (%.1f\u2013%.1f)", 100 * p, 100 * (centre - half), 100 * (centre + half))
}

# The interval alone, "44.2-56.6", for a footnote.
ci_only <- function(x, n) sub("^[^(]*\\((.*)\\)$", "\\1", pct_ci(x, n))

# Percent alone, as a table cell. The tables carry point estimates only; the one interval
# worth reporting, dementia sensitivity (the cell on the fewest people), goes in the
# footnote -- see results_table_for.
pct <- function(x, n) sprintf("%.1f", 100 * x / n)

# `n_tested` / `n_eligible` are supplied by the caller rather than read off `d`, because
# the denominator is row-specific -- see results_table_for.
as_pct_row <- function(d, time, strategy, n_tested, n_eligible) {
  m    <- metric_counts(d)
  cell <- function(k) pct(m[[k]][["x"]], m[[k]][["n"]])
  data.frame(Time                     = time,
             Strategy                 = strategy,
             `Eligible tested, N (%)` = sprintf("%s (%.1f)", format(n_tested, big.mark = ","),
                                                100 * n_tested / n_eligible),
             Sensitivity              = cell("sens"),
             `Sensitivity, MCI`       = cell("sens_mci"),
             `Sensitivity, dementia`  = cell("sens_dem"),
             Specificity              = cell("spec"),
             PPV                      = cell("ppv"),
             NPV                      = cell("npv"),
             Accuracy                 = cell("acc"),
             sens_dem_ci              = ci_only(m$sens_dem[["x"]], m$sens_dem[["n"]]),
             check.names = FALSE)
}

# "First Cycle" is read at each person's own first cycle due for a test, so no calendar
# age is involved. "By Program End" is the calendar snapshot at end_age. Both are
# program-level: eligible people who were never tested count as misses in each.
# anchor = "first_test" swaps in the test-level first row; see anchored_performance for
# why that row is NOT comparable with the program-end one.
#
# "Eligible tested" is the number of eligible people the strategy actually tested, and
# their share. Its
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
               n_tested = fv$n_tested, n_eligible = fv$n_visited)
  }))

  later <- do.call(rbind, lapply(keys, function(scen) {
    st <- strategy_stats[match(scen, strategy_stats$scenario), ]
    as_pct_row(all_test_data %>% filter(scenario == scen, age == at_end_followup),
               "By Program End", reg_row(scen)$label,
               n_tested = st$n_people, n_eligible = st$n_eligible)
  }))

  # The dementia-sensitivity intervals ride along as an attribute rather than a column,
  # for the table's footnote: one row per table row, Time and Strategy to label it.
  out <- rbind(first, later)
  ci  <- out[, c("Time", "Strategy", "sens_dem_ci")]
  out$sens_dem_ci <- NULL
  attr(out, "sens_dem_ci") <- ci
  out
}

# Number needed to test: how many BHAs a strategy runs per 10 people it correctly
# identifies. All four NNT columns keep the SAME numerator -- every test the
# strategy runs -- because no test can be aimed at one stage alone; only the
# denominator changes. So the stage columns read as "tests needed to find 10 people
# at that stage", and they do not average to the overall column.
#
# "Correctly identified" splits three ways, on the person's state at the positive
# test, and the three sum to it:
#   Caught at MCI        impaired with MCI (SEV 0) when they tested positive
#   Caught at dementia   impaired with dementia (SEV >= 1) when they tested positive
#   Caught early         not yet impaired when they tested positive -- outright
#                        healthy or in TCI -- and impaired later in the window.
#                        TCI is SYN 0.5, not impairment, so a TCI-cycle positive
#                        is early. TCI is roughly two thirds of this group in each
#                        main strategy (strategy_stats$n_early_at_tci has the split).
# Early catches used to sit inside "caught at MCI", because they enter the
# identified state at conversion and conversion is by definition MCI. They are
# their own column now so "caught at MCI" means MCI at the time of the test.
#
# The table says "detected" for all of these, and its "correctly detected" is
# narrower than n_identified: only people impaired at the positive test, MCI plus
# dementia. A positive given before onset is reported beside it as a count, with no
# ratio -- the test that gave it was not detecting impairment, which had not begun.
# So the table's overall ratio is NOT strategy_stats$nnt10_any, which counts early
# catches in its denominator.
#
# The shares are another matter: they are of EVERYONE the program flagged who was
# impaired by the end, early catches included, so MCI + dementia + coincidental sum to
# 100 and "any CI" is MCI + dementia. They say what a positive turned out to be; the
# ratios say what it took to find impairment that was already there.
#
# Deliberately not the results table's definition. There a converter is credited to
# program sensitivity from conversion onward, which is accurate: from then on the
# program does hold a correct positive on them. NNT asks something else -- how many
# tests it takes to detect someone while they are impaired -- and an early catch is
# not that. It stays out of the denominator here, and the results table keeps
# crediting it from conversion. Do not align the two.
#
# Ratios are point estimates only, as the paper reports them. Their 95% bootstrap
# intervals stay in strategy_stats (nnt10_*_lo / _hi, from nnt_boot()) if wanted.
#
# Column order is MCI, dementia, then any CI and the early catches: the two stages
# side by side, then the total and what sits outside it. The two timing columns come
# from onset_timing() (see timing_table_for), so they describe a narrower group
# than the counts beside them -- people whose impairment began during the program
# while they were in it. Years from onset is among those caught after onset, which
# leaves out everyone already impaired at the start of the window; years before
# onset is among the early catches, and onset_timing()'s entry rules (a provider,
# no prior clinical diagnosis) drop a handful of them. The Word note says so.
nnt_table_for <- function(keys) {
  if (skip_incomplete(keys, "NNT table")) return(NULL)
  d <- strategy_stats[match(keys, strategy_stats$scenario), ]
  n_detected <- d$n_at_mci + d$n_at_dem
  n_all      <- n_detected + d$n_early_catch   # = n_identified
  one_dp   <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
  timing_mean_sd <- function(outcome, col) vapply(keys, function(scen) {
    t <- scenario_reports[[scen]]$timing
    x <- t[[col]][t$outcome == outcome]
    if (length(x) < 2) return(NA_character_)
    sprintf("%.1f (%.1f)", mean(x), sd(x))
  }, character(1), USE.NAMES = FALSE)
  data.frame(
    Strategy                                           = reg_row(keys)$label,
    `Tests administered`                               = d$n_tests,
    `Tests per person tested`                          = d$per_person,
    `MCI at detection, n (%)`                          = n_pct(d$n_at_mci, n_all),
    `MCI at detection, tests per 10`                    = one_dp(d$nnt10_at_mci),
    `Dementia at detection, n (%)`                     = n_pct(d$n_at_dem, n_all),
    `Dementia at detection, tests per 10`               = one_dp(d$nnt10_at_dem),
    `Any CI at detection, n (%)`                       = n_pct(n_detected, n_all),
    `Any CI at detection, tests per 10`                 = one_dp(d$nnt10_detected),
    `Any CI at detection, years from onset, mean (SD)` = timing_mean_sd("caught", "lag"),
    `Coincidentally detected before onset, n (%)`      = n_pct(d$n_early_catch, n_all),
    `Coincidentally detected before onset, years before onset, mean (SD)` = timing_mean_sd("early", "lead"),
    check.names = FALSE)
}

# Time from impairment onset to identification, from onset_timing(). The four
# outcome columns partition the incident cases. Lag statistics are among those
# caught after onset only; early catches have their own column instead of a lag,
# and their own years-flagged-before-onset, the time they spent labeled while
# unimpaired (1 = flagged the year before onset).
#
# "Within k years" is caught by a positive test no more than k years after the onset
# year, among people whose onset left at least k years of window to be caught in --
# onset at 79 cannot show a 2-year catch by 80, so it would dilute the share. Deaths
# and clinical diagnoses inside those k years count as not caught. Point estimate only,
# as in the results table.
n_pct <- function(x, n) sprintf("%s (%.1f)", format(x, big.mark = ",", trim = TRUE), 100 * x / n)

timing_table_for <- function(keys) {
  if (skip_incomplete(keys, "timing table")) return(NULL)
  do.call(rbind, lapply(keys, function(scen) {
    d   <- scenario_reports[[scen]]$timing
    n   <- nrow(d)
    k_n <- function(o) sum(d$outcome == o)
    lag  <- d$lag[d$outcome == "caught"]
    lead <- d$lead[d$outcome == "early"]
    med_iqr <- function(x) {
      q <- quantile(x, c(0.25, 0.5, 0.75), type = 1, names = FALSE)
      sprintf("%d (%d\u2013%d)", q[2], q[1], q[3])
    }
    mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x), sd(x))
    within <- function(k) {
      at_risk <- d$onset_age <= max(testing_window) - k
      pct(sum(at_risk & d$outcome == "caught" & d$lag <= k, na.rm = TRUE), sum(at_risk))
    }
    data.frame(Strategy                          = reg_row(scen)$label,
               `Impaired during program, N`      = format(n, big.mark = ","),
               `Flagged before onset, n (%)`     = n_pct(k_n("early"), n),
               `Caught after onset, n (%)`       = n_pct(k_n("caught"), n),
               `Clinic diagnosed first, n (%)`   = n_pct(k_n("clinic"), n),
               `Not caught, n (%)`               = n_pct(k_n("neither"), n),
               `Years flagged before onset, median (IQR)` = med_iqr(lead),
               `Years flagged before onset, mean (SD)`    = mean_sd(lead),
               `Years to catch, median (IQR)`    = med_iqr(lag),
               `Years to catch, mean (SD)`       = mean_sd(lag),
               `Caught within 1 year, %`         = within(1),
               `Caught within 2 years, %`        = within(2),
               check.names = FALSE)
  }))
}

timing_table <- timing_table_for(main_keys)
as_flex(timing_table)

# Time spent labeled while unimpaired, from labeled_unimpaired(): the harm side of
# early detection. Person-years run from each person's first positive while
# unimpaired to conversion, death, or end of program, whichever comes first. Split
# by whether the person is ever impaired (at any age, as early_pos is), and scaled
# per 1,000 ever-eligible so strategies of different reach compare.
labeled_table_for <- function(keys) {
  if (skip_incomplete(keys, "labeled table")) return(NULL)
  do.call(rbind, lapply(keys, function(scen) {
    d  <- scenario_reports[[scen]]$labeled
    st <- strategy_stats[match(scen, strategy_stats$scenario), ]
    py <- sum(d$years)
    data.frame(Strategy                            = reg_row(scen)$label,
               `People labeled while unimpaired`   = format(nrow(d), big.mark = ","),
               `Person-years labeled`              = format(py, big.mark = ","),
               `Years per person, mean (SD)`       = sprintf("%.1f (%.1f)", mean(d$years), sd(d$years)),
               `Person-years, later impaired, n (%)` = n_pct(sum(d$years[d$later_impaired]), py),
               `Person-years, never impaired, n (%)` = n_pct(sum(d$years[!d$later_impaired]), py),
               `Person-years per 1,000 eligible`   = format(round(1000 * py / st$n_eligible),
                                                            big.mark = ","),
               check.names = FALSE)
  }))
}

labeled_table <- labeled_table_for(main_keys)
as_flex(labeled_table)

# The paper's table: timing and labeling together, since they answer one question --
# how early does a strategy find people, and what does finding them early cost.
# Transposed, measures down and strategies across, because there are three strategies
# and over a dozen measures. Each duration is given as both median (IQR) and mean (SD):
# years are whole cycles, so the medians bunch at 1-2 and the means carry the spread.
# Block is for the Word report to group rows on.
timing_burden_table_for <- function(keys) {
  tt <- timing_table_for(keys)
  lt <- labeled_table_for(keys)
  if (is.null(tt) || is.null(lt)) return(NULL)
  block <- function(name, src, cols) {
    m <- t(as.matrix(src[, unname(cols), drop = FALSE]))
    colnames(m) <- src$Strategy
    data.frame(Block = name, Measure = names(cols), m, check.names = FALSE, row.names = NULL)
  }
  burden <- "Labeled positive while unimpaired"
  rbind(
    block("Impairment beginning during the program", tt, c(
      "People, N"                         = "Impaired during program, N",
      "Flagged before onset, n (%)"       = "Flagged before onset, n (%)",
      "Detected after onset, n (%)"       = "Caught after onset, n (%)",
      "Clinic diagnosed first, n (%)"     = "Clinic diagnosed first, n (%)",
      "Not detected, n (%)"               = "Not caught, n (%)")),
    block("Time to detection, among those detected after onset", tt, c(
      "Years from onset, median (IQR)"    = "Years to catch, median (IQR)",
      "Years from onset, mean (SD)"       = "Years to catch, mean (SD)",
      "Detected within 1 year, %"  = "Caught within 1 year, %",
      "Detected within 2 years, %" = "Caught within 2 years, %")),
    # Consecutive rows with one Block value make one group, so the early-catch lead
    # time heads the burden block although it comes from the timing table
    block(burden, tt, c(
      "Flagged before onset: years flagged, median (IQR)" = "Years flagged before onset, median (IQR)",
      "Flagged before onset: years flagged, mean (SD)"    = "Years flagged before onset, mean (SD)")),
    block(burden, lt, c(
      "People, N"                         = "People labeled while unimpaired",
      "Person-years"                      = "Person-years labeled",
      "Years per person, mean (SD)"       = "Years per person, mean (SD)",
      "Person-years, later impaired, n (%)" = "Person-years, later impaired, n (%)",
      "Person-years, never impaired, n (%)" = "Person-years, never impaired, n (%)",
      "Person-years per 1,000 eligible"   = "Person-years per 1,000 eligible"))
  )
}

timing_burden_table <- timing_burden_table_for(main_keys_by_intensity)
as_flex(timing_burden_table)

nnt_table <- nnt_table_for(main_keys)
as_flex(nnt_table)

# Same table for the other groups, if the call runs long
nnt_table_sens_selective <- nnt_table_for(sens_selective_keys)
nnt_table_sens_inclusive <- nnt_table_for(sens_inclusive_keys)
nnt_table_pcp            <- nnt_table_for(pcp_keys)


results_table <- results_table_for(main_keys)
as_flex(results_table)
# sensitivity of the program, not sensitivity of the test

results_table_sens_selective <- results_table_for(sens_selective_keys)
as_flex(results_table_sens_selective)

results_table_sens_inclusive <- results_table_for(sens_inclusive_keys)
as_flex(results_table_sens_inclusive)

results_table_pcp <- results_table_for(pcp_keys)
as_flex(results_table_pcp)

