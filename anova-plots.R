unscaled <- test_maker(100, 20, 0.8, 0, 0.4, c(1:10))

library(ggplot2)

# -------------------------------------------------------------------------- #

# ------------------------------------------------------------- #
# Plots mean correlation under each treatment for each variable #
# ------------------------------------------------------------- #

get_m_scores <- function(data, A0 = character(0)) {
  
  A <- A0
  taxa <- rownames(data[[1]])
  n_treat <- length(data)
  
  # row-wise standardized values
  uvw_vals <- lapply(data, function(mat) t(scale(t(mat))))
  
  out <- list()
  
  for (j in taxa) {
    
    Av <- compute_uAvAwA(uvw_vals, A, j)
    
    m_1j <- cov(uvw_vals[[1]][j, ], Av$u_A) # plot the collection of these
    m_2j <- cov(uvw_vals[[2]][j, ], Av$v_A)
    m_3j <- cov(uvw_vals[[3]][j, ], Av$w_A)
    
    out[[j]] <- tibble(
      Taxa = j,
      Treatment = c("T1", "T2", "T3"),
      m_score = c(m_1j, m_2j, m_3j)
    )
  }
  
  bind_rows(out)
}

m_df <- get_m_scores(unscaled, c("V1","V2","V3","V4","V5","V6"))

subset_taxa <- c("V1", "V2", "V3", "V20")

m_df |>
  filter(Taxa %in% subset_taxa) |>
  ggplot(aes(Treatment, m_score, color = Treatment)) +
  geom_point() +
  geom_line(aes(group = Taxa)) +
  facet_wrap(~ Taxa) +
  theme_minimal()

# ---------------------------------------------------------------------------- #

# ------------------------------------------------------------ #
# Plots individual contribution to correlations                #
# Shows skewness for contributions under treatments w high cov #
# ------------------------------------------------------------ #

get_contributions <- function(data, A0 = character(0)) {
  
  A <- A0
  taxa <- rownames(data[[1]])
  
  uvw_vals <- lapply(data, function(mat) t(scale(t(mat))))
  
  out <- list()
  
  for (j in taxa) {
    Av <- compute_uAvAwA(uvw_vals, A, j)
    u_contributions <- uvw_vals[[1]][j, ] * Av$u_A
    v_contributions <- uvw_vals[[2]][j, ] * Av$v_A
    w_contributions <- uvw_vals[[3]][j, ] * Av$w_A
    
    n1 <- length(u_contributions)
    n2 <- length(v_contributions)
    n3 <- length(w_contributions)
    
    out[[j]] <- tibble(
      Taxa = j,
      Treatment = rep(c("T1","T2","T3"), times = c(n1, n2, n3)),
      Value = c(u_contributions, v_contributions, w_contributions)
    )
  }
  bind_rows(out)
}


contributions <- get_contributions(unscaled, 
                      c("V1","V2","V3","V4","V5","V6", "V7", "V8", "V9", "V10"))

subset_taxa <- c("V1", "V2", "V3", "V91")

contributions |>
  filter(Taxa %in% subset_taxa) |>
  ggplot(aes(Treatment, Value, color = Treatment)) +
  geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.8) +
  facet_wrap(~ Taxa) +
  theme_minimal()

# ---------------------------------------------------------------------------- #

# --------------------------------------------- #
# Plots ANOVA-esque situation                   #
# --------------------------------------------- #

contributions |>
  filter(Taxa == "V1") |>
  ggplot(aes(Treatment, Value, color = Treatment)) +
  # individual contributions (u_ij * u_iA)
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  # horizontal mean line (S_k(j,A))
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    fatten = 1,
    linewidth = 0.8
  ) +
  theme_minimal()

contributions |>
  filter(Taxa == "V80") |>
  ggplot(aes(Treatment, Value, color = Treatment)) +
  # individual contributions (u_ij * u_iA)
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  # horizontal mean line (S_k(j,A))
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    fatten = 1,
    linewidth = 0.8
  ) +
  theme_minimal()

# -------------------------------------------- #
# Plots of 60 standardized measurements        #
# note: not that useful, only for confirmation #
# -------------------------------------------- #

make_long_df <- function(data, taxa_subset = NULL) {
  taxa_all <- rownames(data[[1]])
  if (is.null(taxa_subset)) {
    taxa <- taxa_all
  } else {
    taxa <- intersect(taxa_subset, taxa_all)
  }
  long_df <- do.call(
    rbind,
    lapply(seq_along(data), function(k) {
      df <- data[[k]][taxa, , drop = FALSE]
      data.frame(
        taxon = rep(taxa, each = ncol(df)),
        treatment = factor(paste0("T", k)),
        value = as.vector(t(df))
      )
    })
  )
  return(long_df)
}

plot_grid <- function(data, taxa_subset = NULL) {
  df_long <- make_long_df(data, taxa_subset)
  ggplot(df_long, aes(x = treatment, y = value, color = treatment)) +
    geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.8) +
    facet_wrap(~ taxon, scales = "fixed") +
    theme_bw() +
    theme(
      strip.text = element_text(size = 9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = "Treatment",
      y = "Value",
      title = "All 60 measurements per row"
    )
}

plot_grid(unscaled, taxa_subset = c("V1", "V2", "V10", "V20"))
