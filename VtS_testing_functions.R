
compute_S <- function(corr_matrix, j, A) { # where j is the target and A is a set
  A <- setdiff(A, j)
  if (length(A) == 0) return(0)
  return(1/length(A) * (sum(corr_matrix[j,A])))
}

compute_delta <- function(corr_matrices, j, A) {
  s_values <- sapply(corr_matrices, 
                     function(corr_matrix) compute_S(corr_matrix, j, A))
  avg_s <- 1/length(corr_matrices) * sum(s_values)
  return(sum((s_values-avg_s)^2))
}



# Main affinity cycle using compute_delta as the metric
find_A_star_delta <- function(corr_matrices, A0 = character(0), alpha = 0.05, n_perm=1000) {
  taxa <- rownames(corr_matrices[[1]])
  A <- A0
  reps <- 0
  repeat {
    delta_values <- list()
    for (tax in taxa){
      if (tax %in% A) {
        val <- compute_delta(corr_matrices, tax, setdiff(A, tax))
        delta_values <- append(delta_values, val)
      }
      else{
        val <- compute_delta(corr_matrices, tax, A)
        delta_values <- append(delta_values, val)
      }
    }
    names(delta_values) <- taxa
  
  delta_perm <- matrix(NA, nrow = length(delta_values), ncol = n_perm)
  rownames(delta_perm) <- names(delta_values)
  
  # lump data for permutations
  n_groups <- length(corr_matrices)
  combined_data <- do.call(cbind, corr_matrices)
  
  # Create vector of group labels (assumes samples in each original_data_list are columns)
  group_sizes <- sapply(corr_matrices, ncol)
  group_labels <- rep(1:n_groups, times = group_sizes)
  
  for (perm_i in 1:n_perm) {
    # Permute group labels
    permuted_labels <- sample(group_labels)
    
    # Split combined_data columns by permuted groups
    permuted_data_list <- lapply(1:n_groups, function(g) {
      combined_data[, permuted_labels == g, drop = FALSE]
    })
    
    # Compute correlation matrices on permuted data
    corr_matrices_perm <- lapply(permuted_data_list, function(data) cor(t(data)))
    
    # Compute delta per taxa on permuted correlation matrices
    delta_perm[, perm_i] <- sapply(names(delta_values), function(j) {
      compute_delta(corr_matrices_perm, j, A)
    })
  }
  
  # Calculate p-values per taxa (proportion of permuted deltas >= observed delta)
  p_values <- sapply(1:length(delta_values), function(i) {
    mean(delta_perm[i, ] >= delta_values[i])
  })
  names(p_values) <- names(delta_values)
  sig <- names(p_values)[p.adjust(p_values, method = "BH") < alpha]
  if (setequal(A, sig)) break
  print(A)
  A <- sig
  reps = reps + 1
  print(reps)
}
  return(A)
}





# delta_values <- sapply(taxa, function(j) {
#   compute_delta(corr_matrices, j, A)
# })