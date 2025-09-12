library(MASS)
sigma_A <- diag(100)  # makes a 100x100 identity matrix
sigma_A[1:10, 1:10] <- 0.5 #gives the first 10 rows positive correlations
sigma_A[row(sigma_A) == col(sigma_A)] = 1 # makes sure diagonal is one

sample_A <- mvrnorm(100, rep(0, 100), sigma_A)
sample_B <- mvrnorm(100, rep(0, 100), sigma_A)
sample_C <- mvrnorm(100, rep(0, 100), sigma_A)

colnames(sample_A) <- paste0("V", 1:100)
rownames(sample_A) <- paste0("V", 1:100)
colnames(sample_B) <- paste0("V", 1:100)
rownames(sample_B) <- paste0("V", 1:100)
colnames(sample_C) <- paste0("V", 1:100)
rownames(sample_C) <- paste0("V", 1:100)

test_matrices <- lapply(list(sample_A, sample_B, sample_C), cor)

compute_delta(test_matrices, "V1", c("V10"))

compute_delta(test_matrices, "V1", c("V100"))

find_A_star_delta(test_matrices, c("V1"))

find_A_star_delta(test_matrices, c("V4", "V2", "V8"))

