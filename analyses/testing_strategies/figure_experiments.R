######################################## PAPER 2 MAIN FIGURE: DESIGN EXPERIMENTS ########################################
#
# Sandbox for re-designing the main testing-strategies figure. Nothing here is
# sourced by the analysis: it reads the cached per-cycle counts in
# test_perf_results/ and writes candidate figures to plots/experiments/. The
# production plotting code in test_performance_helpers.R is left untouched.
#
#   source("analyses/testing_strategies/figure_experiments.R")
#
# The problem being solved: the published version encodes true cognitive status
# in the LINE colour and test result in the POINT colour, so a false negative
# (red line, blue points) and a false positive (blue line, red points) render as
# near-identical blue-and-red objects wherever they cross. Colour also always
# wins the eye, so whichever series is drawn last reads as the only one there.
# Each variant below moves one of the two facts off colour and onto a channel
# that survives being overdrawn -- dash pattern, marker shape, or a panel of its
# own.

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)


## Inputs ####
run_id   <- "20260825_124029"
perf_dir <- "analyses/testing_strategies/test_perf_results"
out_dir  <- "analyses/testing_strategies/plots/experiments"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_ages <- 65:80

# Least to most intensive, which is also how the paper discusses them
main_keys   <- c("r1bhapos", "s1bhapos_emr", "u3bhapos_rand50")
main_labels <- c(r1bhapos        = "Reactive",
                 s1bhapos_emr    = "Selective",
                 u3bhapos_rand50 = "Inclusive")


## Data ####
# early_pos folds into fp and converted_tp into tp, matching show_early_pos = FALSE
# in the production figure. Deaths are dropped throughout: a cumulative total is
# not a test outcome, and at ~42,000 it sets the y-axis for everything else.
load_counts <- function(keys = main_keys, labels = main_labels,
                        ages = plot_ages, non_testers = FALSE) {
  d <- lapply(keys, function(k) {
    readRDS(file.path(perf_dir, paste0(k, "_", run_id, ".rds"))) %>%
      mutate(scenario = k)
  }) %>%
    bind_rows() %>%
    filter(age %in% ages) %>%
    transmute(Age      = age,
              Scenario = factor(labels[scenario], levels = unname(labels)),
              TP = tp + converted_tp,
              FP = fp + early_pos,
              TN = tn,
              FN = fn,
              `No test, healthy`  = notest_tn,
              `No test, impaired` = notest_fn)

  keep <- if (non_testers) names(d) else setdiff(names(d), c("No test, healthy", "No test, impaired"))
  d %>%
    select(all_of(keep)) %>%
    pivot_longer(-c(Age, Scenario), names_to = "Outcome", values_to = "Count") %>%
    mutate(Outcome = factor(Outcome, levels = c("TP", "FN", "TN", "FP",
                                                "No test, impaired", "No test, healthy")),
           Truth = factor(ifelse(grepl("TP|FN|impaired", Outcome), "Impaired", "Healthy"),
                          levels = c("Impaired", "Healthy")),
           Test  = factor(case_when(Outcome %in% c("TP", "FP") ~ "Positive",
                                    Outcome %in% c("TN", "FN") ~ "Negative",
                                    TRUE                       ~ "Not tested"),
                          levels = c("Positive", "Negative", "Not tested")),
           Agree = factor(ifelse(Outcome %in% c("TP", "TN"), "Test agrees", "Test errs"),
                          levels = c("Test agrees", "Test errs")))
}

d_main    <- load_counts()
d_with_nt <- load_counts(non_testers = TRUE)


## Shared look ####
# Full outcome names, spelled out. "FP" is jargon the reader decodes on every
# glance, and a figure this size has room for words.
outcome_long <- c(TP = "True positive",  FN = "False negative",
                  TN = "True negative",  FP = "False positive",
                  `No test, impaired` = "Not tested, impaired",
                  `No test, healthy`  = "Not tested, healthy")

# Hue carries true status (warm = impaired, cool = healthy); within a hue, the
# dark saturated shade is the test getting it right and the light shade the test
# getting it wrong. The four series separate by hue AND by value, and the two
# error series still read as a pair.
pal_outcome <- c(TP = "#9E2A2B", FN = "#E8A33D",
                 TN = "#1B4965", FP = "#5FA8D3",
                 `No test, impaired` = "#8C6D46", `No test, healthy` = "#9AABB5")

pal_truth <- c(Impaired = "#9E2A2B", Healthy = "#2A6F97")

base_theme <- theme_minimal(base_size = 15) +
  theme(panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        panel.spacing.x    = unit(1.9, "lines"),   # keeps 80 and 65 off each other
        panel.spacing.y    = unit(1.4, "lines"),
        axis.title         = element_text(face = "bold"),
        strip.text         = element_text(size = 15, face = "bold", hjust = 0),
        plot.title         = element_text(face = "bold", size = 17),
        plot.subtitle      = element_text(colour = "grey35", size = 12),
        legend.position    = "bottom",
        legend.title       = element_text(face = "bold"),
        legend.key.width   = unit(1.6, "cm"))

y_count <- function(...) scale_y_continuous(labels = comma, ...)
x_age   <- function(right = 0) {
  scale_x_continuous(breaks = seq(65, 80, 5),
                     expand = expansion(mult = c(0.03, right)))
}
y_sqrt <- function(breaks = c(0, 1000, 5000, 15000, 30000, 50000)) {
  scale_y_sqrt(labels = comma, breaks = breaks)
}

# End-of-line labels, so the reader never travels to a legend and back.
#
# Two things have to be got right or direct labelling is worse than a legend:
# the labels must all land in ONE panel (facets repeat the data, so an ungrouped
# slice_max silently draws them in every panel at once), and labels for series
# that converge must be pushed apart. `spread_labels` walks the series from the
# bottom up and forces a minimum vertical gap, working in the transformed space
# so it behaves on a square-root axis too.
spread_labels <- function(y, gap, trans = c("identity", "sqrt")) {
  trans <- match.arg(trans)
  f  <- if (trans == "sqrt") sqrt else identity
  fi <- if (trans == "sqrt") function(x) x^2 else identity
  ord <- order(y)
  z <- f(y[ord])

  # Resolve collisions by splitting the difference between the two labels
  # involved, then relaxing until nothing moves. A single bottom-up pass instead
  # anchors the lowest label and shunts everything above it upward, which is why
  # labels used to float above the line ends they belong to.
  for (iter in seq_len(80)) {
    moved <- FALSE
    for (i in seq_along(z)[-1]) {
      deficit <- gap - (z[i] - z[i - 1])
      if (deficit > 1e-9) {
        z[i - 1] <- z[i - 1] - deficit / 2
        z[i]     <- z[i]     + deficit / 2
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  out <- numeric(length(y))
  out[ord] <- fi(pmax(z, 0))
  out
}

end_labels <- function(d, row_var = NULL, size = 4, gap_frac = 0.06,
                       trans = "identity", panel = c("last", "first", "all"), nudge = 0.4) {
  panel <- match.arg(panel)

  # "all" labels every panel; otherwise pick one so the labels are not repeated
  # three times over. Labelling every panel needs the same right-hand expansion
  # in each, which is why the x scale below is shared.
  if (panel != "all") {
    scen <- levels(d$Scenario)[if (panel == "last") nlevels(d$Scenario) else 1]
    d    <- d %>% filter(Scenario == scen)
  }

  # A single shared y-axis is just the one-group case, so give it a constant
  # grouping column rather than a NULL `by`, which would join on nothing.
  d  <- d %>% mutate(.row = "all")
  gv <- if (is.null(row_var)) ".row" else row_var

  # free_y in facet_grid is shared down a row, so the gap is a fraction of that
  # row's full range across all scenarios -- not of the labelled panel alone.
  rng <- d %>%
    group_by(across(all_of(gv))) %>%   # free_y is shared down a row, not per panel
    summarise(span = {
      f <- if (trans == "sqrt") sqrt else identity
      f(max(Count)) - f(min(Count))
    }, .groups = "drop")

  lab <- d %>%
    filter(Age == max(Age)) %>%
    left_join(rng, by = gv) %>%
    group_by(across(all_of(c("Scenario", gv)))) %>%
    mutate(y = spread_labels(Count, gap_frac * first(span), trans),
           lab = outcome_long[as.character(Outcome)]) %>%
    ungroup()

  geom_text(data = lab, aes(x = Age, y = y, label = lab), hjust = 0,
            nudge_x = nudge, size = size, fontface = "bold",
            show.legend = FALSE, inherit.aes = FALSE,
            colour = pal_outcome[as.character(lab$Outcome)])
}


## V0 -- published encoding, non-tester lines removed ####
# The literal ask: the figure as it stands, minus the two "No test" series that
# were setting the y-axis. Kept as the before/after reference -- the colour
# collision it was drawn with is still here.
v0 <- ggplot(d_main, aes(Age, Count, group = Outcome)) +
  geom_line(aes(colour = Truth), linewidth = 0.9) +
  geom_point(aes(fill = Test), shape = 21, colour = "white", size = 2.4, stroke = 0.6) +
  facet_wrap(~Scenario) +
  scale_colour_manual("Cognitive status", values = pal_truth) +
  scale_fill_manual("Test result", values = c(Positive = "#9E2A2B", Negative = "#2A6F97")) +
  x_age() + y_count() +
  labs(title = "V0. Published encoding, non-tester lines removed",
       subtitle = "Line colour = true status, point colour = test result. Unchanged otherwise",
       y = "People") +
  base_theme

## V1 -- one colour per outcome, labelled at the line end ####
# Truth and test result stop competing for the colour channel: each series is
# one thing, with one colour and its own name at the end of the line.
v1 <- ggplot(d_main, aes(Age, Count, colour = Outcome, group = Outcome)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.7) +
  end_labels(d_main, size = 4) +
  facet_wrap(~Scenario) +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  x_age(right = 0.5) + y_count() +
  labs(title = "V1. One colour per outcome, labelled directly",
       subtitle = "Warm = impaired, cool = healthy; dark = the test was right, light = the test was wrong",
       y = "People") +
  base_theme

## V2 -- hue = truth, dash + shade = whether the test was right ####
# Colour keeps its original job (who the person actually is). Being wrong is a
# second, independent channel: dashed and lighter. Dash pattern is still legible
# through an overdrawn line, which is exactly where the old version failed.
v2 <- ggplot(d_main, aes(Age, Count, colour = Outcome, linetype = Agree, group = Outcome)) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~Scenario) +
  scale_colour_manual("Outcome", values = pal_outcome,
                      breaks = c("TP", "FN", "TN", "FP"),
                      labels = unname(outcome_long[c("TP", "FN", "TN", "FP")])) +
  scale_linetype_manual("Test verdict", values = c("Test agrees" = "solid", "Test errs" = "22")) +
  x_age() + y_count() +
  labs(title = "V2. Hue = true status; dashed and lighter = the test got it wrong",
       subtitle = "No point markers at all -- overlap is resolved by dash pattern",
       y = "People") +
  base_theme +
  guides(colour = guide_legend(order = 1, nrow = 2,
                               override.aes = list(linetype = "solid", linewidth = 1.4)),
         linetype = guide_legend(order = 2, nrow = 2,
                                 override.aes = list(colour = "grey25", linewidth = 1)))

## V3 -- hue = truth, marker shape = test result ####
# Points take their own line's colour, so nothing is bicoloured and there is no
# discordance to misread. The test result lives entirely in the marker: solid
# disc = positive, hollow disc = negative. Shape reads correctly whichever
# series happens to be on top.
v3 <- ggplot(d_main, aes(Age, Count, colour = Truth, group = Outcome)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(aes(shape = Test), size = 2.9, stroke = 1.1, fill = "white") +
  facet_wrap(~Scenario) +
  scale_colour_manual("Cognitive status", values = pal_truth) +
  scale_shape_manual("Test result", values = c(Positive = 16, Negative = 21)) +
  x_age() + y_count() +
  labs(title = "V3. Hue = true status, marker shape = test result",
       subtitle = "Filled = tested positive, hollow = tested negative. A point never disagrees with its line",
       y = "People") +
  base_theme +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4)))

## V4 -- split the panel by true cognitive status ####
# The structural fix. A false negative and a false positive cannot overlap
# because they are drawn in different rows. Free y per row lets the impaired
# row -- an order of magnitude smaller -- be read at all.
v4 <- ggplot(d_main, aes(Age, Count, colour = Outcome, group = Outcome)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.7) +
  end_labels(d_main, row_var = "Truth", size = 3.9) +
  facet_grid(Truth ~ Scenario, scales = "free_y", switch = "y") +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  x_age(right = 0.62) + y_count(n.breaks = 5) +
  labs(title = "V4. One row per true cognitive status, free y-axis",
       subtitle = "Top row: people who really are impaired. Bottom row: people who are healthy. The rows are on different scales",
       y = "People") +
  base_theme +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

## V5 -- four colours on a shared square-root y-axis ####
# One scale for all three strategies, so cross-panel comparison stays honest,
# but compressed at the top so the reactive panel is not a flat line on the
# floor. The single biggest legibility win of anything here.
v5 <- ggplot(d_main, aes(Age, Count, colour = Outcome, group = Outcome)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.7) +
  end_labels(d_main, size = 3.9, trans = "sqrt") +
  facet_wrap(~Scenario) +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  x_age(right = 0.52) + y_sqrt() +
  labs(title = "V5. One colour per outcome, shared square-root y-axis",
       subtitle = "One scale across strategies; the root transform keeps the small series off the floor",
       y = "People (square-root scale)") +
  base_theme

## V6 -- composition of the tested population ####
# A different question -- the MIX of outcomes rather than the counts -- and
# immune to overlap by construction. Good as a companion panel, not a
# replacement: it hides how many people each strategy tests at all.
v6 <- ggplot(d_main, aes(Age, Count, fill = Outcome)) +
  geom_area(position = "fill", colour = "white", linewidth = 0.25) +
  facet_wrap(~Scenario) +
  scale_fill_manual("Outcome", values = pal_outcome, labels = outcome_long) +
  x_age() + scale_y_continuous(labels = percent) +
  labs(title = "V6. Share of everyone tested, by outcome",
       subtitle = "Composition rather than counts. Nothing can overlap, but the size of each programme is lost",
       y = "Share of tests") +
  base_theme +
  guides(fill = guide_legend(nrow = 1))

## V7 -- errors separated from correct calls ####
# Correct calls are the bulk; the errors are the finding. Two rows on their own
# scales stops the large series from flattening the small ones.
d_split <- d_main %>%
  mutate(Panel = factor(ifelse(Agree == "Test agrees", "Test got it right", "Test got it wrong"),
                        levels = c("Test got it right", "Test got it wrong")))
v7 <- ggplot(d_split, aes(Age, Count, colour = Outcome, group = Outcome)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.7) +
  end_labels(d_split, row_var = "Panel", size = 3.9) +
  facet_grid(Panel ~ Scenario, scales = "free_y", switch = "y") +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  x_age(right = 0.62) + y_count(n.breaks = 5) +
  labs(title = "V7. Correct calls and errors in separate rows, free y-axis",
       subtitle = "The error row is the finding; on a shared axis it is squashed against zero",
       y = "People") +
  base_theme +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 90))

## V8 -- V5's encoding with the non-tester series added back ####
# For checking whether the redesign survives if the untested have to stay in.
# They are an order of magnitude larger, so the root axis is doing real work.
v8 <- ggplot(d_with_nt, aes(Age, Count, colour = Outcome, group = Outcome)) +
  geom_line(aes(linetype = Test == "Not tested"), linewidth = 1.05) +
  end_labels(d_with_nt, size = 3.6, trans = "sqrt") +
  facet_wrap(~Scenario) +
  scale_colour_manual(values = pal_outcome, guide = "none") +
  scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "21"), guide = "none") +
  x_age(right = 0.62) + y_sqrt(c(0, 1000, 5000, 15000, 30000, 60000)) +
  labs(title = "V8. Non-tested series kept, square-root axis",
       subtitle = "Dotted = never tested. For comparison with the versions that drop them",
       y = "People (square-root scale)") +
  base_theme

## V9 -- shape encoding on the square-root axis ####
# V3's marker idea plus V5's scale: the two changes are independent, and the
# root axis is what finally separates the small series in the reactive panel.
v9 <- ggplot(d_main, aes(Age, Count, colour = Truth, group = Outcome)) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(aes(shape = Test), size = 2.9, stroke = 1.1, fill = "white") +
  facet_wrap(~Scenario) +
  scale_colour_manual("Cognitive status", values = pal_truth) +
  scale_shape_manual("Test result", values = c(Positive = 16, Negative = 21)) +
  x_age() + y_sqrt() +
  labs(title = "V9. Marker shape for the test result, square-root y-axis",
       subtitle = "V3's encoding on V5's scale",
       y = "People (square-root scale)") +
  base_theme +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.4)))

## V10 -- dash encoding on the square-root axis ####
# V2's dash idea plus V5's scale, for the same reason.
v10 <- ggplot(d_main, aes(Age, Count, colour = Outcome, linetype = Agree, group = Outcome)) +
  geom_line(linewidth = 1.15) +
  facet_wrap(~Scenario) +
  scale_colour_manual("Outcome", values = pal_outcome,
                      breaks = c("TP", "FN", "TN", "FP"),
                      labels = unname(outcome_long[c("TP", "FN", "TN", "FP")])) +
  scale_linetype_manual("Test verdict", values = c("Test agrees" = "solid", "Test errs" = "22")) +
  x_age() + y_sqrt() +
  labs(title = "V10. Dashed = the test got it wrong, square-root y-axis",
       subtitle = "V2's encoding on V5's scale",
       y = "People (square-root scale)") +
  base_theme +
  guides(colour = guide_legend(order = 1, nrow = 2,
                               override.aes = list(linetype = "solid", linewidth = 1.4)),
         linetype = guide_legend(order = 2, nrow = 2,
                                 override.aes = list(colour = "grey25", linewidth = 1)))


## V11 -- V8 over V6, counts above composition ####
# Two rows answering the two halves of the question: how many people land in
# each cell (top, counts, square-root axis so the small series survive), and
# what share of the tested population each cell holds (bottom, composition).
# Every series is labelled at the line end in every panel, so neither row needs
# a legend.

# Programme size per strategy, read from the raw arrays and cached since that is
# the only thing here that needs them. n_people is the denominator of the bottom
# row (all four series are cummax stocks, so someone tested at 65 is still in it
# at 80); n_tests is test volume, which is the larger and faster-growing number.
n_tested_file <- file.path(out_dir, "n_tested.rds")
if (file.exists(n_tested_file)) {
  n_stats <- readRDS(n_tested_file)
} else {
  cycles <- plot_ages - 50 + 1
  n_stats <- do.call(rbind, lapply(main_keys, function(k) {
    a <- readRDS(sprintf("analyses/testing_strategies/sim_results/scenario_%s_sim_%s.rds", k, run_id))$output
    bha <- a[cycles, "BHA", ]
    out <- data.frame(scenario = k,
                      n_people = sum(apply(bha >= 0, 2, any, na.rm = TRUE)),
                      n_tests  = sum(bha >= 0, na.rm = TRUE))
    rm(a, bha); invisible(gc())
    out
  }))
  n_stats$per_person <- n_stats$n_tests / n_stats$n_people
  saveRDS(n_stats, n_tested_file)
}

strip_stats <- setNames(
  sprintf("%s ever tested   \u00b7   %s tests   \u00b7   %.1f per person",
          comma(n_stats$n_people), comma(n_stats$n_tests), n_stats$per_person),
  n_stats$Scenario)

# Both rows share this expansion so the two panel grids line up and 65-80 sits
# at the same horizontal position in each.
x_shared <- x_age(right = 0.46)

## V11 now lives in the production helper, so this sandbox does not carry a
# second copy to drift out of step with it. Tune the figure there.
source("analyses/testing_strategies/test_performance_helpers.R")

v11 <- plot_counts_and_shares(
  do.call(rbind, lapply(main_keys, function(k) {
    readRDS(file.path(perf_dir, paste0(k, "_", run_id, ".rds"))) %>% mutate(scenario = k)
  })),
  scenario_names = main_labels[main_keys],
  strategy_stats = n_stats,
  ages           = plot_ages)




## V12 / V13 -- composition as a share of the COHORT, not of those tested ####
#
# PI's question: instead of normalising the bottom row to 100% of the tested, let
# the stack total be the share of the cohort tested, so growth in reach is visible
# and the untested / deceased remainder shows as area rather than being divided out.
#
# Closing the stack to 100% needs one band the perf data does not carry. At any
# age every member of the starting cohort is exactly one of:
#
#   dead                                             -> death
#   alive, ever tested                               -> TP / FN / TN / FP
#   alive, never tested, ELIGIBLE                    -> notest_tn / notest_fn
#   alive, never tested, NOT eligible                -> the missing band
#
# That last group is 4,000-7,000 people, split roughly evenly between having a
# clinical diagnosis and having no provider. It is left as its own band rather
# than folded into "not tested", because notest_tn / notest_fn are deliberately
# eligibility-gated -- they measure a coverage gap the programme could still
# close, and someone ineligible is not that. Verified to sum to exactly 100,000.
cohort_file <- file.path(out_dir, "cohort_bands.rds")
if (file.exists(cohort_file)) {
  cohort_bands <- readRDS(cohort_file)
} else {
  cycles <- plot_ages - 50 + 1
  cohort_bands <- do.call(rbind, lapply(main_keys, function(k) {
    a <- readRDS(sprintf("analyses/testing_strategies/sim_results/scenario_%s_sim_%s.rds", k, run_id))$output
    bha <- a[,"BHA",]; dx <- a[,"DX",]; hc <- a[,"HCARE",]; alive <- a[,"ALIVE",]
    et <- apply(bha, 2, cummax) >= 0
    dxlag <- rbind(rep(NA, ncol(dx)), dx[-nrow(dx), , drop = FALSE])
    inelig <- (alive == 1) & !et & !(dxlag == 0 & hc == 1)
    out <- data.frame(scenario = k, age = plot_ages,
                      not_eligible = rowSums(inelig[cycles, ], na.rm = TRUE))
    rm(a); invisible(gc())
    out
  }))
  saveRDS(cohort_bands, cohort_file)
}

pal_cohort <- c(pal_outcome[c("TP", "FN", "TN", "FP")],
                `Not tested`          = "#8FA8A2",
                `Not tested, healthy` = "#6F9A94",
                `Not tested, impaired`= "#8C6D46",
                `Not eligible`        = "#B5AEA4",
                Deceased              = "#D8DCDE")

# split_untested = FALSE gives one "Not tested" band; TRUE splits it healthy /
# impaired, which is the same cut the counts row already uses.
plot_cohort_shares <- function(split_untested = FALSE, label_size = 3.0) {
  wide <- lapply(main_keys, function(k) {
    readRDS(file.path(perf_dir, paste0(k, "_", run_id, ".rds"))) %>% mutate(scenario = k)
  }) %>% bind_rows() %>%
    filter(age %in% plot_ages) %>%
    left_join(cohort_bands, by = c("scenario", "age")) %>%
    transmute(Age = age,
              Scenario = factor(main_labels[scenario], levels = unname(main_labels[main_keys])),
              TP = tp + converted_tp, FP = fp + early_pos, TN = tn, FN = fn,
              `Not tested, healthy`  = notest_tn,
              `Not tested, impaired` = notest_fn,
              `Not eligible`         = not_eligible,
              Deceased               = death)

  if (!split_untested) {
    wide <- wide %>%
      mutate(`Not tested` = `Not tested, healthy` + `Not tested, impaired`,
             .keep = "unused")
  }

  # position_stack puts the FIRST factor level on TOP, so the level order is the
  # reverse of the bottom-up reading order.
  untested_lv <- if (split_untested) c("Not tested, healthy", "Not tested, impaired") else "Not tested"
  lv <- c("Deceased", "Not eligible", rev(untested_lv), "FP", "TN", "FN", "TP")

  d <- wide %>%
    pivot_longer(-c(Age, Scenario), names_to = "Band", values_to = "Count") %>%
    mutate(Band = factor(Band, levels = lv)) %>%
    group_by(Scenario, Age) %>%
    mutate(Share = Count / sum(Count)) %>%     # sums to exactly 1: verified partition
    ungroup()

  lab <- d %>%
    filter(Age == max(Age)) %>%
    group_by(Scenario) %>%
    arrange(desc(Band), .by_group = TRUE) %>%
    mutate(ymid = cumsum(Share) - Share / 2,
           y    = f.spread_labels(ymid, 0.062),
           lab  = ifelse(as.character(Band) %in% names(outcome_long),
                         outcome_long[as.character(Band)], as.character(Band))) %>%
    ungroup()

  ggplot(d, aes(Age, Share, fill = Band)) +
    geom_area(colour = "white", linewidth = 0.25) +
    geom_text(data = lab, aes(x = Age, y = y, label = lab), hjust = 0, nudge_x = 0.3,
              size = label_size, fontface = "bold", inherit.aes = FALSE,
              colour = pal_cohort[as.character(lab$Band)]) +
    facet_wrap(~Scenario) +
    scale_fill_manual(values = pal_cohort, guide = "none") +
    x_age(right = 0.62) +
    scale_y_continuous(labels = percent, expand = expansion(mult = c(0.035, 0.02))) +
    labs(y = "Share of the age-50 cohort") +
    base_theme
}

v12 <- plot_cohort_shares(split_untested = FALSE) +
  labs(title = "V12. Share of the whole cohort, untested as one band",
       subtitle = "Stack totals 100% of the starting cohort. The tested block at the bottom is programme reach")

v13 <- plot_cohort_shares(split_untested = TRUE) +
  labs(title = "V13. Share of the whole cohort, untested split by true status",
       subtitle = "As V12, but the eligible-untested band is cut healthy vs impaired")


## Render ####
variants <- list(v0 = v0, v1 = v1, v2 = v2, v3 = v3, v4 = v4, v5 = v5,
                 v6 = v6, v7 = v7, v8 = v8, v9 = v9, v10 = v10, v11 = v11, v12 = v12, v13 = v13)

sizes <- list(v0 = c(13, 6),   v1 = c(14.5, 6),  v2 = c(13, 6.8), v3 = c(13, 6.5),
              v4 = c(15, 8.5), v5 = c(14.5, 6),  v6 = c(13, 6.3), v7 = c(15, 8.5),
              v8 = c(15, 6.5), v9 = c(13, 6.5),  v10 = c(13, 6.8),
              v11 = c(17, 12.5), v12 = c(15, 7), v13 = c(15, 7))

for (nm in names(variants)) {
  ggsave(file.path(out_dir, paste0(nm, ".png")), variants[[nm]],
         width = sizes[[nm]][1], height = sizes[[nm]][2], dpi = 160, bg = "white")
}

cat("Wrote", length(variants), "variants to", out_dir, "\n")
