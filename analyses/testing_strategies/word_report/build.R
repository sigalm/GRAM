######################################## BUILD THE WORD REPORT ########################################
# Run from the project root:  Rscript analyses/testing_strategies/word_report/build.R
#
# Steps, in order:
#   00  reference_landscape.docx   page setup pandoc knits against
#   01  report_data.rds + counts_and_shares.png   paper 2 tables and main figure, rebuilt
#       from an existing simulation run (nothing is re-simulated)
#   Rmd assembles all of it into GRAM_testing_strategies_figures_and_tables.docx
#
# 01 is the slow step; comment it out when its outputs are already current. The
# calibration/validation figures (now docs/_s2_figures.R) are no longer part of this
# report, which covers the testing strategies only.
#
# Paths are written out in full at each call rather than held in a variable: step 01
# sources model/setup.R, which opens with rm(list = ls()), so nothing assigned in
# this file survives them.

# RStudio ships pandoc but does not put it on PATH for a plain Rscript session.
if (!nzchar(Sys.getenv("RSTUDIO_PANDOC"))) {
  local({
    candidates <- file.path("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
                            c("aarch64", "x86_64"))
    found <- candidates[file.exists(file.path(candidates, "pandoc"))]
    if (length(found)) Sys.setenv(RSTUDIO_PANDOC = found[1])
  })
}

source("analyses/testing_strategies/word_report/00_reference_docx.R")
source("analyses/testing_strategies/word_report/01_paper2_tables.R")

rmarkdown::render("analyses/testing_strategies/word_report/gram_testing_strategies_figures_and_tables.Rmd",
                  output_file = "GRAM_testing_strategies_figures_and_tables.docx",
                  quiet = TRUE)

cat("\nWrote analyses/testing_strategies/word_report/GRAM_testing_strategies_figures_and_tables.docx\n")
