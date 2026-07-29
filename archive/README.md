# Archive

Files kept for the record but no longer part of the working model. Nothing here is sourced,
read, or referenced by live code — check before reviving anything, since the surrounding code
has moved on.

| File | Why it is here |
|------|----------------|
| `benchmarking_v3.Rmd` | The pre-protocol calibration notebook. Superseded by `calibration/run_calibration.R`, which implements the formal protocol. Unrunnable as written: it sources `calibration/gram_benchmarking_helpers.R` and reads `data/acs_data/acs_sample_1.rds` / `acs_sample_2.rds`, none of which exist. Kept whole rather than repaired, because its value is the narrative record of which parameter adjustments were tried and why. |
| `allcause_mortality_bysex_clean.RDS` | Intermediate all-cause mortality rates by sex and age. The `saveRDS()` that produced it was removed when `data/mortality/calculate_non_dementia_death_rates.Rmd` was rewritten to carry sex long and pivot once at the end, so the file no longer has a generator. Nothing reads it. |
