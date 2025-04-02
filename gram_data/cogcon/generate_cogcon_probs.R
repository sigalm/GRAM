## Probaility of subjective cognitive concerns by age and cognitive status

m.cogcon_spon <- data.frame(age = 50:100,
                       h = c(0, rep(0.01, 10), rep(0.03, 10), rep(0.05, 10), rep(0.07, 10), rep(0.10, 10)),
                       mci = rep(0.6, 51),
                       dem = rep(0.9, 51))

m.cogcon_elic <- m.cogcon_spon %>%
  mutate(h = h + 0.1,
         mci = mci + 0.1,
         dem = dem + 0.1)


saveRDS(m.cogcon_spon, "gram_data/cogcon/m.cogcon_spon.RDS")
saveRDS(m.cogcon_elic, "gram_data/cogcon/m.cogcon_elic.RDS")




