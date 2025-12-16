## Distribution of medical burden by age
# Multimorbid condition count -- continuous variable 0-15
# Assume exponential


rate <-  0.05
age <- 50:99
prob <- 1 - exp(-rate * (age-50))
hist(prob)

prob <- pmin(pmax(prob, 0), 1)


random <- runif(1000)
medbur <- qbinom(random, size = 15, prob = prob)

plot(age, prob)
hist(medbur)
plot(age,medbur)


# Starting MEDBUR - these values achieve approximately a mean of 2 and variance of 1.5
samples <- round(rbeta(n=100, 0.4, 18) * 15, 0)
hist(samples)
prop.table(table(samples))
mean(samples)

samples2 <- round(rbeta(n = 100, shape1 = 2, shape2 = 18) * 15,0)
hist(samples2)
prop.table(table(samples2))
mean(samples2)


# Starting EDU
mean_desired <- 10.9
sd_desired <- 5
max_value <- 25

# Scaled mean and variance
scaled_mean <- mean_desired / max_value
scaled_variance <- (sd_desired^2) / (max_value^2)

# Calculate alpha and beta
alpha <- scaled_mean * ((scaled_mean * (1 - scaled_mean)) / scaled_variance - 1)
beta <- (1 - scaled_mean) * ((scaled_mean * (1 - scaled_mean)) / scaled_variance - 1)


samples <- rbeta(100, alpha, beta) * max_value
hist(samples)
mean(samples)
sd(samples)


