# Generate question endorsement rates by true cognitive status, P(endorses | state)
# Source: GRAMish workbook, "P conditional testing" tab, 14 Sep 2026.
#
# The workbook gives P(tested | state) for the question-based arm: endorsing a concern AND
# taking the test up, combined. It was built from three observed facts (15% of patients
# endorse a concern, 30% of those take the test up, 55% of those tested are impaired after
# correcting for misclassification), an assumed dementia:MCI testing ratio of 4 at that 30%
# uptake, and a scale-up to the 50% base case, where the ratio falls to 2.9:
#
#   P(tested | state) at 50% uptake   0.041 healthy / 0.342 MCI / 0.986 dementia
#
# The model splits testing into selection (probs_select) and take-up (probs_accept). The
# question arm takes the eRADAR arm's take-up unchanged, so endorsement is the workbook's
# P(tested) divided by the eRADAR take-up at the same 50% base case. Endorsing is then a
# fixed property of the question, as the flag is of eRADAR, and outreach acts only through
# take-up.
#
# Split at 50% rather than 30%: the workbook's 30% figures (0.03 / 0.225 / 0.90) test 5.4%
# of GRAM's population against the 4.5% observed, while the 50% figures test 7.46% against
# 7.5%. Checked back against the observed facts at the eRADAR 30% take-up, these rates give
# 14.3% endorsing, 32.3% of them tested, 56.0% impaired among the tested (observed 15 / 30 /
# 55). Dementia take-up is 1 at every uptake, so P(tested | dementia) is the endorsement
# rate itself at both levels: 0.986 at 30% too, where the workbook had 0.90.
#
# Depends on the eRADAR take-up: regenerate this file whenever the eRADAR flag file or the
# eRADAR config's uptake / p_impaired_tested change. The prevalences are the eRADAR ones,
# saved alongside for the same reason they are saved there.

source("model/helpers/program_utility_helpers.R")

p_tested_50 <- c(h = 0.041, mci = 0.342, dem = 0.986)   # workbook, 50% uptake

eradar <- readRDS("data/eradar/eradar_flag_by_state.RDS")
eradar_accept <- f.accept_matrix(f.select_matrix(h = eradar$h, mci = eradar$mci, dem = eradar$dem),
                                 uptake            = 0.50,   # eRADAR config base case
                                 p_impaired_tested = 0.45,   # eRADAR config base case
                                 prev_mci          = eradar$prev_mci,
                                 prev_dem          = eradar$prev_dem)[1, c("h", "mci", "dem")]

endorse <- p_tested_50 / unlist(eradar_accept)
stopifnot(all(endorse >= 0), all(endorse <= 1))

question_endorse <- data.frame(h = endorse[["h"]], mci = endorse[["mci"]], dem = endorse[["dem"]],
                               prev_mci = eradar$prev_mci, prev_dem = eradar$prev_dem)

print(round(question_endorse, 4))   # 0.1099 / 0.4130 / 0.986

saveRDS(question_endorse, "data/question/question_endorse_by_state.RDS")
