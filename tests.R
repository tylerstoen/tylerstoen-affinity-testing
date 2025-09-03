sigma_A <- diag(100)  # makes a 100x100 identity matrix
sigma_A[1:10, 1:10] <- 0.5 #gives the first 10 rows positive correlations
sigma_A[row(sigma_A) == col(sigma_A)] = 1 # makes sure diagonal is one

sample_A <- mvrnorm(100, mu = 0, sigma = sigma_A)
sample_B <- mvrnorm(100)
sample_C <- mvrnorm(100)
