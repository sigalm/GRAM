# Generate eRADAR flag rates by true cognitive status, P(flagged | state)
# From the published eRADAR validation at the 85th-percentile cutoff (citation to be added).
#
#   dem   the published dementia sensitivity.
#   mci   backed out of the sensitivity with MCI and dementia both counted as cases:
#         sens_impaired = share_dem x sens_dem + (1 - share_dem) x P(flag | MCI), where
#         share_dem is dementia's share of those cases.
#   h     whatever fills the rest of the flagged fraction. The cutoff is a program
#         decision, so the top 15% stay flagged, and specificity follows from that rather
#         than being taken from the paper (87.8% here; 84.7% published).
#
# The paper's reference standard is a dementia diagnosis within the next 12 months; GRAM's
# flag acts on true current status. The two differ in both directions: the paper's cases
# include people who convert during the year (MCI in GRAM when scored, presumably harder
# to flag), and GRAM's dementia includes people who never come to clinical attention
# (presumably less visible in the record). With no clear net direction, the published
# rates are taken as they are. They are the instrument's ceiling: no level of uptake can
# test more of a state than the flag selects.
#
# The translation needs the population the cutoff is applied in. These are GRAM's
# undiagnosed prevalences at age 70 before the July 2026 recalibration. The eRADAR sample's
# mean age is 73, where the current model gives 7.63% MCI / 1.03% dementia; those give
# near-identical rates, so these are kept. They are saved alongside the rates because the
# eRADAR acceptance (f.accept_matrix() in the selective configs) has to be built in the
# same population.

sens_dem      <- 0.659          # P(flag | dementia)
sens_impaired <- 0.577          # P(flag | MCI or dementia)
share_dem     <- 7631 / 11255   # dementia diagnoses among MCI + dementia diagnoses
flag_share    <- 0.15           # the 85th-percentile cutoff
prev_mci      <- 0.0825
prev_dem      <- 0.0093

sens_mci <- (sens_impaired - share_dem * sens_dem) / (1 - share_dem)
p_h      <- (flag_share - sens_mci * prev_mci - sens_dem * prev_dem) / (1 - prev_mci - prev_dem)
stopifnot(sens_mci >= 0, sens_mci <= 1, p_h >= 0, p_h <= 1)

eradar_flag <- data.frame(h = p_h, mci = sens_mci, dem = sens_dem,
                          prev_mci = prev_mci, prev_dem = prev_dem)

print(round(eradar_flag, 4))   # 0.1217 / 0.4043 / 0.659
cat("Implied specificity (MCI and dementia as cases):", round(1 - p_h, 3), "\n")

saveRDS(eradar_flag, "data/eradar/eradar_flag_by_state.RDS")
