# calibration/write_calibrated_params.R
#
# Promotes a calibration run's best-fit point into model/config/calibrated_params.R,
# the file model/setup.R reads to calibrate every run by default.
#
# Kept separate from run_calibration.R so it can be sourced on its own and pointed at
# an already-saved results file, without re-running the several-hour grid search.
#
# The search grid is the archive; the generated file is the promoted point estimate.
# Promotion is deliberately an explicit step rather than setup.R reaching into the
# grid itself: that keeps the values a plain-text, reviewable diff, and stops a
# re-run of the calibration from silently changing the model under an analysis.

f.write_calibrated_params <- function(results,
                                      out_file = "model/config/calibrated_params.R",
                                      results_file = NA_character_,
                                      digits = 4) {

  best <- results[!results$error, ]
  best <- best[which.min(best$gof), ]

  prov <- attr(results, "provenance")
  if (is.null(prov)) prov <- list()   # results saved before provenance was added

  # Grid geometry is recovered from the results themselves rather than taken as an
  # argument, so this works the same on a live run and on a file read back from disk.
  axis <- function(p) sort(unique(results[[p]]))
  rng  <- function(p) range(axis(p))
  step <- function(p) { a <- axis(p); if (length(a) < 2) NA_real_ else mean(diff(a)) }

  fmt_num <- function(x) formatC(x, format = "f", digits = digits)
  fmt_vec <- function(x) sprintf("c(%s)", paste(fmt_num(x), collapse = ", "))
  fmt_chr <- function(x) if (is.null(x) || is.na(x)) "NA" else sprintf('"%s"', x)

  lines <- c(
    '# ****************************************************************',
    '# ======= GRAM CALIBRATED PARAMETERS =======',
    '# ****************************************************************',
    '#',
    '# GENERATED FILE -- DO NOT EDIT BY HAND.',
    '# Written by calibration/write_calibrated_params.R.',
    '# To change these values, re-run the calibration and let it rewrite this file.',
    '#',
    '# This is the single source of truth for GRAM\'s calibrated parameters. model/setup.R',
    '# sources it and merges l.calibrated into l.inputs, so every run is calibrated by',
    '# default and nothing needs to be applied at the call site.',
    '#',
    sprintf('# Values are rounded to %d decimal places, which is far finer than the grid', digits),
    '# resolution the search actually resolves (see param_step below), so the rounding',
    '# discards no information. The unrounded best-fit row is kept in best_unrounded.',
    '',
    'l.calibrated <- list(',
    sprintf('  param1  = %s,   # multiplier on the published age-specific MCI incidence hazards',
            fmt_num(best$param1)),
    sprintf('  param2a = %s,   # curvature of the CDR-SB progression age curve during MCI',
            fmt_num(best$param2a)),
    sprintf('  param2b = %s    # scales the maximum CDR-SB progression rate during MCI',
            fmt_num(best$param2b)),
    ')',
    '',
    'attr(l.calibrated, "provenance") <- list(',
    sprintf('  run_id         = %s,', fmt_chr(prov$run_id)),
    sprintf('  calibrated_at  = %s,', fmt_chr(prov$saved_at)),
    sprintf('  git_commit     = %s,', fmt_chr(prov$git_commit)),
    sprintf('  results_file   = %s,', fmt_chr(results_file)),
    sprintf('  written_at     = %s,', fmt_chr(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    sprintf('  n_calib        = %s,', if (is.null(prov$n_calib)) "NA" else format(prov$n_calib)),
    sprintf('  n_grid         = %d,', nrow(results)),
    sprintf('  gof            = %s,', fmt_num(best$gof)),
    sprintf('  gof_components = list(wssd_mci = %s, wssd_dem = %s, wssd_mort = %s,',
            fmt_num(best$wssd_mci), fmt_num(best$wssd_dem), fmt_num(best$wssd_mort)),
    sprintf('                        n_mci = %d, n_dem = %d, n_mort = %d),',
            best$n_mci, best$n_dem, best$n_mort),
    sprintf('  param_range    = list(param1  = %s,', fmt_vec(rng("param1"))),
    sprintf('                        param2a = %s,', fmt_vec(rng("param2a"))),
    sprintf('                        param2b = %s),', fmt_vec(rng("param2b"))),
    sprintf('  param_step     = list(param1 = %s, param2a = %s, param2b = %s),',
            fmt_num(step("param1")), fmt_num(step("param2a")), fmt_num(step("param2b"))),
    sprintf('  best_unrounded = list(param1 = %.10g, param2a = %.10g, param2b = %.10g)',
            best$param1, best$param2a, best$param2b),
    ')',
    ''
  )

  writeLines(lines, out_file)

  cat(sprintf("\n=== WROTE: %s ===\n", out_file))
  cat(sprintf("  param1  = %s\n", fmt_num(best$param1)))
  cat(sprintf("  param2a = %s\n", fmt_num(best$param2a)))
  cat(sprintf("  param2b = %s\n", fmt_num(best$param2b)))
  cat(sprintf("  from run_id %s (commit %s), gof = %s\n",
              prov$run_id, prov$git_commit, fmt_num(best$gof)))

  invisible(best)
}
