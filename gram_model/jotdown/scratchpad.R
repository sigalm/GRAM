

# Change in cognition (MMSE) based on PET-based Braak staging
# Ref: Biel 2021 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8361801/)


mmse_change_B0 <- rnorm(1, mean=-0.08, sd=0.27)
mmse_change_B1 <- rnorm(1, mean=-0.39, sd=0.56)
mmse_change_B4 <- rnorm(1, mean=-0.59, sd=0.51)
mmse_change_B6 <- rnorm(1, mean=-1.60, sd=1.01)
mmse_change_Batypical <- rnorm(1, mean=-0.55, sd=0.80)



# Cross-sectional association in functional scores (CDR-SB) and PET-based Braak stage
# Note SE = estimate / t-value
# Ref: Macedo 2024 

cdr_sb_f <- b0 + b1*Braak12 + b2*Braak34 + b3*Braak56 + bage*age + bsex*sex + bsuvr*suvr

b0 <- rnorm(1, mean=0.05, sd=0.05/0.39)
b1 <- rnorm(1, mean=0.15, sd=0.15/0.89)
b2 <- rnorm(1, mean=1.17, sd=1.17/6.18)
bage <- rnorm(1, mean=-0.09, sd=0.09/1.80)
bsex <- rnorm(1, mean=0.07, sd=0.07/0.71)
bsuvr <- rnorm(1, mean=0.16, sd=0.16/2.26)

