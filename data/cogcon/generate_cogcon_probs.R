## Probability of subjective cognitive concerns by age and cognitive status

# From KP prelim data:
# 2% of population gets referred to geriatrics / BHA (assume equivalent to spontaneous)
# 15% without prior dx gets flagged for BHA in annual welness visit
# We do not have breakdown by actual diagnosis, but assume more likely referral with increased severity



m.cogcon_spon <- data.frame(age = 50:100,
                            h = rep(0.01, 51),
                            mci = rep(0.10, 51),
                            dem = rep(0.30, 51))

m.cogcon_elic <- data.frame(age = 50:100,
                            h = rep(0.05, 51),
                            mci = rep(0.40, 51),
                            dem = rep(0.80, 51))


saveRDS(m.cogcon_spon, "gram_data/cogcon/m.cogcon_spon.RDS")
saveRDS(m.cogcon_elic, "gram_data/cogcon/m.cogcon_elic.RDS")




