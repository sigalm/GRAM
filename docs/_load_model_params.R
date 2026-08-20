# docs/_load_model_params.R
#
# Lets the supplements quote model parameters by reading them from the model, so a
# parameter change cannot leave the prose describing a version that no longer exists.
# (docs/01 spent a calibration cycle quoting param2a = 1.80 / param2b = 1.668 and a
# progression curve that had stopped existing; that is what this is here to prevent.)
#
# Usage, in a chunk near the top of the document:
#
#   ```{r model-params, include=FALSE}
#   source(if (file.exists("docs/_load_model_params.R")) "docs/_load_model_params.R"
#          else "_load_model_params.R")
#   inp  <- gram_inputs()
#   prov <- gram_provenance()
#   ```
#
# then inline, e.g.:  `r fmt(inp$param2a)`
#
# Underscore-prefixed so rmarkdown::render_site() and knitr treat it as a helper
# rather than a document to render.

# Resolve the repository root whether the document is knitted from docs/ (knitr's
# default working directory is the document's own folder) or from the project root.
gram_root <- function() {
  if (file.exists("model/setup.R")) return(normalizePath("."))
  if (file.exists("../model/setup.R")) return(normalizePath(".."))
  stop("Cannot locate the GRAM repository root from ", getwd())
}

# Returns the full l.inputs list, calibrated, exactly as a simulation would see it.
#
# setup.R is sourced into a throwaway environment rather than the document's own. It
# opens with rm(list = ls()), which at the top level of a sourced file resolves to the
# environment it is being evaluated in -- so with sys.source(envir = ) that clears the
# throwaway environment and leaves the knit environment alone. Sourcing it normally
# would wipe the document mid-knit.
gram_inputs <- function() {
  env <- new.env()
  old <- setwd(gram_root())
  on.exit(setwd(old), add = TRUE)
  suppressMessages(sys.source("model/setup.R", envir = env))
  env$l.inputs
}

# Provenance of the calibration run the current values came from: run_id, git commit,
# date, GOF, searched range. Reads the generated params file directly, which is cheap
# and dependency-free, so it can be called without loading the whole model.
gram_provenance <- function() {
  env <- new.env()
  old <- setwd(gram_root())
  on.exit(setwd(old), add = TRUE)
  sys.source("model/config/calibrated_params.R", envir = env)
  attr(env$l.calibrated, "provenance")
}

# The calibrated CDR-SB progression age curve for slow progressors, ages 50-100.
# Derived with the model's own f.CDRslow_curve() definition where available; the
# fallback keeps this file usable without sourcing the helper files.
gram_CDRslow_curve <- function(inp) {
  (seq(0, 1, length.out = 51)^inp$param2a) * (inp$param2b * inp$r.CDRslow_mean)
}

# Fixed-decimal formatter for inline use; avoids scientific notation and stray
# trailing digits in prose.
fmt <- function(x, digits = 2) formatC(x, format = "f", digits = digits)

# Prints which calibration the document was built from, for whoever is knitting it.
#
# These supplements are public-facing, so run ids and git commits do not belong in the
# rendered output -- but they are exactly what you want when checking internally that a
# doc was built against the calibration you think it was. stderr is the one channel that
# achieves both: knitr captures a chunk's stdout (discarded by results='hide') and its
# messages (discarded by message=FALSE, both implied by include=FALSE), but writes stderr
# straight through to the console. Verified: nothing below reaches the knitted file.
gram_announce_provenance <- function(prov = gram_provenance()) {
  cat(sep = "", file = stderr(),
      "\n--- GRAM calibration provenance (console only; not in the rendered document) ---\n",
      "  run_id      : ", prov$run_id,        "\n",
      "  git commit  : ", prov$git_commit,    "\n",
      "  calibrated  : ", prov$calibrated_at, "\n",
      "  GOF         : ", prov$gof,           "\n",
      "  results file: ", prov$results_file,  "\n",
      "-------------------------------------------------------------------------------\n\n")
  invisible(prov)
}
