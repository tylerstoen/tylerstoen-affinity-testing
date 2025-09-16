library(MASS)

# Test 1 (First 10 correlated differently across treatments, others indp)
test1_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test1_sigma_A[1:10, 1:10] <- 0.5 #gives the first 10 rows positive correlations
test1_sigma_A[row(test1_sigma_A) == col(test1_sigma_A)] = 1 # makes sure diagonal is one

test1_sigma_B <- diag(100)
test1_sigma_B[1:10,1:10] <- 0.8
test1_sigma_B[row(test1_sigma_B) == col(test1_sigma_B)] = 1

test1_sigma_C <- diag(100)
test1_sigma_C[1:10,1:10] <- 0.2
test1_sigma_C[row(test1_sigma_C) == col(test1_sigma_C)] = 1

test1_sample_A <- mvrnorm(100, rep(0, 100), test1_sigma_A)
test1_sample_B <- mvrnorm(100, rep(0, 100), test1_sigma_B)
test1_sample_C <- mvrnorm(100, rep(0, 100), test1_sigma_C)

colnames(test1_sample_A) <- paste0("V", 1:100)
rownames(test1_sample_A) <- paste0("V", 1:100)
colnames(test1_sample_B) <- paste0("V", 1:100)
rownames(test1_sample_B) <- paste0("V", 1:100)
colnames(test1_sample_C) <- paste0("V", 1:100)
rownames(test1_sample_C) <- paste0("V", 1:100)

test1_matrices <- lapply(list(test1_sample_A, test1_sample_B, test1_sample_C), cor)

compute_delta(test1_matrices, "V1", c("V10"))
compute_delta(test1_matrices, "V1", c("V100"))
find_A_star_delta(test1_matrices, c("V1"))
find_A_star_delta(test1_matrices, c("V4", "V2", "V8"))


# Test 2 (First 10 only correlated in sample A)
test2_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test2_sigma_A[1:10, 1:10] <- 0.7 #gives the first 10 rows positive correlations (NOTE 0.5 signal not strong enough)
test2_sigma_A[row(test2_sigma_A) == col(test2_sigma_A)] = 1 # makes sure diagonal is one

test2_sample_A <- mvrnorm(100, rep(0, 100), test2_sigma_A)
test2_sample_B <- mvrnorm(100, rep(0, 100), diag(100))
test2_sample_C <- mvrnorm(100, rep(0, 100), diag(100))

colnames(test2_sample_A) <- paste0("V", 1:100)
rownames(test2_sample_A) <- paste0("V", 1:100)
colnames(test2_sample_B) <- paste0("V", 1:100)
rownames(test2_sample_B) <- paste0("V", 1:100)
colnames(test2_sample_C) <- paste0("V", 1:100)
rownames(test2_sample_C) <- paste0("V", 1:100)

test2_matrices <- lapply(list(test2_sample_A, test2_sample_B, test2_sample_C), cor)

compute_delta(test2_matrices, "V1", c("V8"))
compute_delta(test2_matrices, "V1", c("V18"))
compute_delta(test2_matrices, "V25", c("V80"))
find_A_star_delta(test2_matrices, c("V1")) # Fails at 0.5, works at 0.9
find_A_star_delta(test2_matrices, c("V1", "V2", "V8"))


# Test 3 (Correlated under two conditions, not under the third)
test3_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test3_sigma_A[1:10, 1:10] <- 0.7 #gives the first 10 rows positive correlations (NOTE 0.5 signal not strong enough)
test3_sigma_A[row(test3_sigma_A) == col(test3_sigma_A)] = 1 # makes sure diagonal is one

test3_sigma_B <- diag(100)
test3_sigma_B[1:10, 1:10] <- 0.7
test3_sigma_B[row(test3_sigma_B) == col(test3_sigma_B)] = 1

test3_sample_A <- mvrnorm(100, rep(0, 100), test3_sigma_A)
test3_sample_B <- mvrnorm(100, rep(0, 100), test3_sigma_B)
test3_sample_C <- mvrnorm(100, rep(0, 100), diag(100))

colnames(test3_sample_A) <- paste0("V", 1:100)
rownames(test3_sample_A) <- paste0("V", 1:100)
colnames(test3_sample_B) <- paste0("V", 1:100)
rownames(test3_sample_B) <- paste0("V", 1:100)
colnames(test3_sample_C) <- paste0("V", 1:100)
rownames(test3_sample_C) <- paste0("V", 1:100)

test3_matrices <- lapply(list(test3_sample_A, test3_sample_B, test3_sample_C), cor)

compute_delta(test3_matrices, "V1", c("V8"))
compute_delta(test3_matrices, "V1", c("V18"))
compute_delta(test3_matrices, "V25", c("V80"))
find_A_star_delta(test3_matrices, c("V1"))
find_A_star_delta(test3_matrices, c("V1", "V2", "V8"))
find_A_star_delta(test3_matrices, c("V11", "V20"))


# Test 4
test4_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test4_sigma_A[1:10, 1:10] <- 0.7 #gives the first 10 rows positive correlations (NOTE 0.5 signal not strong enough)
test4_sigma_A[row(test3_sigma_A) == col(test3_sigma_A)] = 1 # makes sure diagonal is one

test4_sigma_B <- diag(100)
test4_sigma_B[11:20, 11:20] <- 0.7
test4_sigma_B[row(test4_sigma_B) == col(test4_sigma_B)] = 1

test4_sample_A <- mvrnorm(100, rep(0, 100), test4_sigma_A)
test4_sample_B <- mvrnorm(100, rep(0, 100), test4_sigma_B)
test4_sample_C <- mvrnorm(100, rep(0, 100), diag(100))

colnames(test4_sample_A) <- paste0("V", 1:100)
rownames(test4_sample_A) <- paste0("V", 1:100)
colnames(test4_sample_B) <- paste0("V", 1:100)
rownames(test4_sample_B) <- paste0("V", 1:100)
colnames(test4_sample_C) <- paste0("V", 1:100)
rownames(test4_sample_C) <- paste0("V", 1:100)

test4_matrices <- lapply(list(test4_sample_A, test4_sample_B, test4_sample_C), cor)

compute_delta(test4_matrices, "V1", c("V8"))
compute_delta(test4_matrices, "V1", c("V18"))
compute_delta(test4_matrices, "V25", c("V80"))
find_A_star_delta(test4_matrices, c("V1"))
find_A_star_delta(test4_matrices, c("V1", "V2", "V8"))
find_A_star_delta(test4_matrices, c("V11", "V20"))
find_A_star_delta(test4_matrices, c("V1", "V11")) # returns all V1-V20
find_A_star_delta(test4_matrices, c("V1", "V2", "V3", "V12"))
find_A_star_delta(test4_matrices, c("V1", "V3", "V4", "V5", "V8", "V9", "V16","V19"))
find_A_star_delta(test4_matrices, c("V80", "V43", "V29"))

# Test 5 (Weaker Signal)
test5_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test5_sigma_A[1:10, 1:10] <- 0.2 #gives the first 10 rows positive correlations (NOTE 0.5 signal not strong enough)
test5_sigma_A[row(test5_sigma_A) == col(test5_sigma_A)] = 1 # makes sure diagonal is one

test5_sigma_B <- diag(100)
test5_sigma_B[1:10, 1:10] <- 0.2
test5_sigma_B[row(test5_sigma_B) == col(test5_sigma_B)] = 1

test5_sample_A <- mvrnorm(100, rep(0, 100), test5_sigma_A)
test5_sample_B <- mvrnorm(100, rep(0, 100), test5_sigma_B)
test5_sample_C <- mvrnorm(100, rep(0, 100), diag(100))

colnames(test5_sample_A) <- paste0("V", 1:100)
rownames(test5_sample_A) <- paste0("V", 1:100)
colnames(test5_sample_B) <- paste0("V", 1:100)
rownames(test5_sample_B) <- paste0("V", 1:100)
colnames(test5_sample_C) <- paste0("V", 1:100)
rownames(test5_sample_C) <- paste0("V", 1:100)

test5_matrices <- lapply(list(test5_sample_A, test5_sample_B, test5_sample_C), cor)

compute_delta(test5_matrices, "V1", c("V8"))
compute_delta(test5_matrices, "V1", c("V18"))
compute_delta(test5_matrices, "V25", c("V80"))
find_A_star_delta(test5_matrices, c("V1"))
find_A_star_delta(test5_matrices, c("V1", "V2", "V8"))


# Test 6 (high corr under all conditions)
test6_sigma_A <- diag(100)  # makes a 100x100 identity matrix
test6_sigma_A[1:10, 1:10] <- 0.7 #gives the first 10 rows positive correlations (NOTE 0.5 signal not strong enough)
test6_sigma_A[row(test6_sigma_A) == col(test6_sigma_A)] = 1 # makes sure diagonal is one

test6_sigma_B <- diag(100)
test6_sigma_B[1:10, 1:10] <- 0.7
test6_sigma_B[row(test6_sigma_B) == col(test6_sigma_B)] = 1

test6_sigma_C <- diag(100)
test6_sigma_C[1:10, 1:10] <- 0.7
test6_sigma_C[row(test6_sigma_C) == col(test6_sigma_C)] = 1

test6_sample_A <- mvrnorm(100, rep(0, 100), test6_sigma_A)
test6_sample_B <- mvrnorm(100, rep(0, 100), test6_sigma_B)
test6_sample_C <- mvrnorm(100, rep(0, 100), test6_sigma_C)

colnames(test6_sample_A) <- paste0("V", 1:100)
rownames(test6_sample_A) <- paste0("V", 1:100)
colnames(test6_sample_B) <- paste0("V", 1:100)
rownames(test6_sample_B) <- paste0("V", 1:100)
colnames(test6_sample_C) <- paste0("V", 1:100)
rownames(test6_sample_C) <- paste0("V", 1:100)

test6_matrices <- lapply(list(test6_sample_A, test6_sample_B, test6_sample_C), cor)

compute_delta(test6_matrices, "V1", c("V8"))
compute_delta(test6_matrices, "V1", c("V18"))
compute_delta(test6_matrices, "V25", c("V80"))
find_A_star_delta(test6_matrices, c("V1"))
find_A_star_delta(test6_matrices, c("V1", "V2", "V8"))