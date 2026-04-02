

find_A_star_delta2 <- function(data, A0 = character(0), alpha = 0.05, n_perm=1000) {
  
  permed_corr_mats_list <- permute_data(data, n_perm)
  
  # cleaned_mat <- vineyard_data_to_matrix(data)
  
  #corr_matrices <- lapply(data, function(mat) cor(t(mat)))
  
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
    
    # filtered_data <- lapply(data, function(df){df %>% filter(rownames(df) %in% taxa)})
      
    for (perm in 1:length(permed_corr_mats_list)) {
      delta_perm[, perm] <- sapply(names(delta_values), function(j) {
        compute_delta(permed_corr_mats_list[[perm]], j, A)
      })
    }
    
    # Calculate p-values per taxa (proportion of permuted deltas >= observed delta)
    p_values <- sapply(1:length(delta_values), function(i) {
      mean(delta_perm[i, ] >= delta_values[i])
    })
    names(p_values) <- names(delta_values)
    p.adj <- p.adjust(p_values, method = "BH")
    sig <- names(p_values)[p.adj < alpha]
    if (setequal(A, sig)) break
    
    
    # Stopping condition for bouncing results
    sig_str <- paste(sort(sig), collapse = ",")
    history_strings <- sapply(A_history, paste, collapse = ",")
    
    if (sig_str %in% history_strings) {
      message("Exiting early: detected a repeating pattern in A.")
      break
    } else{
      A_history[[length(A_history) + 1]] <- sort(sig)
    }
    
    print(A)
    A <- sig
    reps = reps + 1
    print(reps)
    
  }
  return(A)
}


get_permutations <- function(data) {
  n_conditions <- length(data)
  n_row <- nrow(data[[1]])
  n_samples <- ncol(data[[1]][-c(1,2)])
  classes <- data[[1]]$Class
  sample_names <- colnames(data[[1]])[-c(1,2)]
  
  # build permuted dataframes
  permuted_data <- vector("list", n_conditions)
  for (i in seq_len(n_conditions)) {
    permuted_data[[i]] <- data[[i]]
  }
    
    # loop over variables and samples
  for (row in seq_len(n_row)) {
    for (sample_j in seq_len(n_samples)) {
      # extract counts across conditions for this (var, sample)
      counts <- sapply(data, function(df) df[row, sample_names[sample_j]])
        
      # permute across conditions
      permuted_counts <- sample(counts, length(counts))
        
      # save into permuted dataframes
      for (cond_k in seq_len(n_conditions)) {
        permuted_data[[cond_k]][row, sample_names[sample_j]] <- permuted_counts[cond_k]
      }
    }
  }
  return(permuted_data)
}


permute_data <- function(data, nperm) {
  permed_corr_mats <- vector("list", nperm)
  
  for (n in 1:nperm) {
    orig_order <- c(rep("A", ncol(data[[1]])), rep("B", ncol(data[[2]])), rep("C", ncol(data[[3]])))
    all_data_together <- cbind(data[[1]], data[[2]], data[[3]])
    
    perm_1 <- sample(orig_order)
    
    fake_A <- all_data_together[, perm_1 == "A"]
    fake_B <- all_data_together[, perm_1 == "B"]
    fake_C <- all_data_together[, perm_1 == "C"]
    
    rownames(fake_A) <- rownames(data[[1]])
    rownames(fake_B) <- rownames(data[[2]])
    rownames(fake_C) <- rownames(data[[3]])
    
    permed_data <- list(fake_A, fake_B, fake_C)
    permed_corr_mats[[n]] <- lapply(permed_data, function(mat) cor(t(mat)))
  }
  return(permed_corr_mats)
} # returns a list of three correlation matrices





vineyard_data_to_matrix <- function(data) { #input list of data (object from split() works)
  names(data) <- seq_along(data)
  
  # Three dataframe inputs w/ variable, group, and all samples as columns
  mats <- lapply(data, function(df){
    df %>%
      column_to_rownames(var = colnames(df)[1]) %>%  # first column as rownames
      .[, -1]
  }) # drop variable and group cols
  
  zero_var_rows <- unique(unlist(lapply(mats, function(mat) {
    rownames(mat)[apply(mat, 1, var) == 0]
  })))
  
  clean_mats <- lapply(mats, function(mat) {
    mat[!(rownames(mat) %in% zero_var_rows), , drop = FALSE]
  })
  
  return(clean_mats)
}

