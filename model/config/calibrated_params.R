# ****************************************************************
# ======= GRAM CALIBRATED PARAMETERS =======
# ****************************************************************
#
# GENERATED FILE -- DO NOT EDIT BY HAND.
# Written by calibration/write_calibrated_params.R.
# To change these values, re-run the calibration and let it rewrite this file.
#
# This is the single source of truth for GRAM's calibrated parameters. model/setup.R
# sources it and merges l.calibrated into l.inputs, so every run is calibrated by
# default and nothing needs to be applied at the call site.
#
# Values are rounded to 4 decimal places, which is far finer than the grid
# resolution the search actually resolves (see param_step below), so the rounding
# discards no information. The unrounded best-fit row is kept in best_unrounded.

l.calibrated <- list(
  param1  = 1.6408,   # multiplier on the published age-specific MCI incidence hazards
  param2a = 1.9500,   # curvature of the CDR-SB progression age curve during MCI
  param2b = 1.8997    # scales the maximum CDR-SB progression rate during MCI
)

attr(l.calibrated, "provenance") <- list(
  run_id         = "20260814_162215",
  calibrated_at  = "2026-08-14 21:40:49",
  git_commit     = "faa8f32",
  results_file   = "calibration/calibration_results_20260814_162215_pass2.RDS",
  written_at     = "2026-08-20 14:35:08",
  n_calib        = 10000,
  n_grid         = 1000,
  gof            = 21.3122,
  gof_components = list(wssd_mci = 5.2201, wssd_dem = 18.3433, wssd_mort = 825.7495,
                        n_mci = 3, n_dem = 6, n_mort = 50),
  param_range    = list(param1  = c(1.3425, 2.2375),
                        param2a = c(1.3500, 2.2500),
                        param2b = c(1.2510, 2.0850)),
  param_step     = list(param1 = 0.0994, param2a = 0.1000, param2b = 0.0927),
  best_unrounded = list(param1 = 1.640833333, param2a = 1.95, param2b = 1.899666667)
)

