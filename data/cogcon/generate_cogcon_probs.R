## Probability of subjective cognitive concerns by age and cognitive status

# From KP prelim data:
# 2% of population gets referred to geriatrics / BHA (assume equivalent to reactive)
# 15% without prior dx gets flagged for BHA in annual wellness visit (selective)
# We do not have breakdown by actual diagnosis, but assume more likely referral with increased severity

# Reactive: patient-initiated (e.g., self-referral due to concerns)
# Selective: provider-initiated (e.g., annual wellness visit screening)

m.cogcon_reactive <- data.frame(age = 50:100,
                                h = rep(0.01, 51),
                                mci = rep(0.10, 51),
                                dem = rep(0.30, 51))

m.cogcon_selective <- data.frame(age = 50:100,
                                 h = rep(0.05, 51),
                                 mci = rep(0.40, 51),
                                 dem = rep(0.80, 51))


saveRDS(m.cogcon_reactive, "data/cogcon/m.cogcon_reactive.RDS")
saveRDS(m.cogcon_selective, "data/cogcon/m.cogcon_selective.RDS")




