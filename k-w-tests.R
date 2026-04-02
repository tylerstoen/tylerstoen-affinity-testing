# K-W tests

# compute matrix

A_star_k.w <- function(data,
                       A0 = character(0),
                       alpha = 0.05) {
  
  taxa <- rownames(data[[1]])
  
  n1 <- length(data[[1]])
  n2 <- length(data[[2]])
  n3 <- length(data[[3]])
  
  uvw <- compute_uvw(data)
  
  # ----------------------------------
  # NEW INITIALIZATION
  # ----------------------------------
  
  if (length(A0) == 0) {
    A <- initialize_A(uvw, taxa, n1, n2, n3)
  } else {
    A <- A0
  }
  
  reps <- 0
  A_history <- list()
  
  repeat {
    
    pvals <- compute_kw_pvals(uvw, A, taxa, n1, n2, n3)
    
    p.adj <- p.adjust(pvals, method = "BH")
    sig <- names(pvals)[!is.na(p.adj) & p.adj < alpha]
    
    if (length(A) == 1) {
      sig <- union(sig, A)
    }
    
    if (setequal(A, sig)) break
    
    # -----------------------------
    # Bouncing detection
    # -----------------------------
    
    sig_str <- paste(sort(sig), collapse = ",")
    history_strings <- sapply(A_history, paste, collapse = ",")
    
    if (sig_str %in% history_strings) {
      message("Exiting early: detected a repeating pattern in A.")
      break
    } else {
      A_history[[length(A_history) + 1]] <- sort(sig)
    }
    
    print(A)
    A <- sig
    reps <- reps + 1
    print(reps)
  }
  
  return(A)
}

compute_kw_pvals <- function(uvw, A, taxa, n1, n2, n3) {
  
  pvals <- numeric(length(taxa))
  names(pvals) <- taxa
  
  for (j in taxa) {
    
    Av <- compute_uAvAwA(uvw, A, j)
    
    u_tild <- uvw[[1]][j,] * Av$u_A
    v_tild <- uvw[[2]][j,] * Av$v_A
    w_tild <- uvw[[3]][j,] * Av$w_A
    
    x <- c(u_tild, v_tild, w_tild)
    g <- c(rep(1, n1), rep(2, n2), rep(3, n3))
    
    pvals[j] <- kruskal.test(x = x, g = g)$p.value
  }
  
  return(pvals)
}


initialize_A_1 <- function(uvw, taxa, n1, n2, n3, K=10) { # single point of largest difference
  scores <- sapply(taxa, function(j) {
    means <- sapply(uvw, function(mat) mean(mat[j, ]))
    var(means)   # or max difference
  })
  
  start <- names(which.max(scores))
  
  # Step 2: expand using your method
  pvals <- compute_kw_pvals(uvw, A = start, taxa, n1, n2, n3)
  
  # Step 3: take top K
  A_init <- names(sort(pvals))[1:K]
  
  return(A_init)
}

initialize_A <- function(uvw, taxa, n1, n2, n3,
                         K = 10,
                         M = 20) {
  
  scores <- sapply(taxa, function(j) {
    means <- sapply(uvw, function(mat) mean(mat[j, ]))
    var(means)
  })
  
  top_candidates <- names(sort(scores, decreasing = TRUE))[1:min(M, length(taxa))]
  
  best_A <- NULL
  best_score <- Inf
  
  for (start in top_candidates) {
    
    pvals <- compute_kw_pvals(uvw, A = start, taxa, n1, n2, n3)
    
    A_candidate <- names(sort(pvals))[1:K]

    candidate_pvals <- compute_kw_pvals(uvw, A = A_candidate, taxa, n1, n2, n3)
    
    score <- mean(sort(candidate_pvals)[1:K])
    
    if (score < best_score) {
      best_score <- score
      best_A <- A_candidate
    }
  }
  
  return(best_A)
}

 unscaled <- test_maker(100, 20, 0.8, 0, 0.1, c(1:10))

A_star_k.w(unscaled)



