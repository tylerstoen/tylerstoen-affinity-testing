# ---------------------------------------------------------
# FAILED
# A_star <- function(data, A0 = character(0)) {
#   
#   corr_matrices <- lapply(data, function(mat) cor(t(mat)))
#   A <- A0
#   n_treat <- length(corr_matrices)
#   nA <- length(A)
#   
#   F_values <- numeric(nrow(corr_matrices[[1]]))
#   
#   for (j in 1:nrow(corr_matrices[[1]])) {
#     S_vals <- sapply(seq_len(n_treat), function(k) {
#       compute_S(corr_matrices[[k]], j, setdiff(A, j))
#     })
#     
#     S_bar <- mean(S_vals)
#     between <- sum((S_vals - S_bar)^2) * (nA / (n_treat - 1))
#     
#     WSS <- 0
#     for (k in seq_len(n_treat)) {
#       corr_mat <- corr_matrices[[k]]
#       for (i in A) {
#         WSS <- WSS + (corr_mat[j, i] - S_vals[k])^2
#       }
#     }
#     within <- WSS / (nA * n_treat - n_treat)
#     
#     F_values[j] <- between / within
#   }
#   
#   return(F_values)
# }


# -------------------------------------------------------------------------
A_star_uvw <- function(data, A0 = character(0)) {
  
  A <- A0
  taxa <- rownames(data[[1]])
  n_treat <- length(data)
  
  F_vals <- list()
  
  # row-wise standardized values
  uvw_vals <- compute_uvw(data)
  
  for (j in taxa) {
    
    # get u_A, v_A, w_A from helper
    Av <- compute_uAvAwA(uvw_vals, A, j) # these are m's in the write up
    
    u_A <- Av$u_A # m1j
    v_A <- Av$v_A # m2j
    w_A <- Av$w_A # m3j
    
    # projection / alignment scores
    m_1j <- cov(uvw_vals[[1]][j, ], u_A) # S1(i, A)s in the ANOVA write up
    m_2j <- cov(uvw_vals[[2]][j, ], v_A) # S2(i, A)
    m_3j <- cov(uvw_vals[[3]][j, ], w_A) # S3(i, A)
    
    # within sum of squares
    WSS <- sum(
      (uvw_vals[[1]][j, ] * u_A - m_1j)^2 +
        (uvw_vals[[2]][j, ] * v_A - m_2j)^2 +
        (uvw_vals[[3]][j, ] * w_A - m_3j)^2
    )
    
    # between sum of squares
    m_dot_j <- mean(c(m_1j, m_2j, m_3j))
    
    BSS <- sum(
      (m_1j - m_dot_j)^2,
      (m_2j - m_dot_j)^2,
      (m_3j - m_dot_j)^2
    )
    
    n <- sum(vapply(data, ncol, integer(1)))
    
    within <- WSS / (n - n_treat)
    between <- BSS / (n_treat - 1)
    
    F_vals[[j]] <- between / within
  }
  
  p_vals <- lapply(F_vals, function(val) {
    1 - pf(val, n_treat - 1, n - n_treat)
  })
  
  p_vals
}

# flip axes, scale and flip back


compute_uAvAwA <- function(uvw_vals, A, j) {
  
  n_treat <- length(uvw_vals)
  
  # helper: compute column means safely
  col_mean_or_zero <- function(mat, rows) {
    if (length(rows) == 0) {
      rep(0, ncol(mat))
    } else {
      colMeans(mat[rows, , drop = FALSE])
    }
  }
  
  if (j %in% A) {
    A_use <- setdiff(A, j)
  } else {
    A_use <- A
  }
  
  u_A <- col_mean_or_zero(uvw_vals[[1]], A_use)
  v_A <- col_mean_or_zero(uvw_vals[[2]], A_use)
  w_A <- col_mean_or_zero(uvw_vals[[3]], A_use)
  
  list(
    u_A = u_A,
    v_A = v_A,
    w_A = w_A
  )
}


compute_uvw <- function(data) {
  uvw_vals <- lapply(data, function(mat){
    t(scale(t(mat)))
  })
}






unscaled <- test_maker(100, 20, 0.7, 0, 0, c(1:10))


test <- A_star_uvw(unscaled, c("V1"))
