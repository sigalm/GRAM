# Calculate SD from 95% CI to generate random distribution

# Data from Jamalian et al, Figure 3. Using the CREAD2 cohort data
# CDR-SB score change from baseline ranges from 0 - 2.5 over 18 months (mean)
# 5th percentile == -0.1 - 2.2
# 95th percentile == 0.5 - 3.5

range_mean <- (2.5 - 0) * 12/18
range_min <- (2.2 - (-0.1)) * 12/18
range_max <- (3.5 - 0.5) * 12/18 

n <- 398 # see table 2

# sd <- (upper_ci - lower_ci)*sqrt(n)/1.96

sd <- (range_max - range_min)*sqrt(n)/1.96

# mean 1.67, range 4.75 points (not percent!!)

rnorm(n=5, mean = 1.67, sd = 4.75)



