## tests 2
library(tidyverse)
library(MASS)
library(mvtnorm)

## Vineyard test (Doesn't work because of missing values or zero-variance in the data.)
class_dat <- read.csv(here::here("rarified_counts_by_Class.csv"))

class_dat <- class_dat %>% 
  mutate(
    sample = paste(Year, Season, Replicate)
  )

cc_pivoted <- class_dat %>% 
  filter(Experiment == "Crop Cover") %>% 
  group_by(Class, sample, Treatment) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = total,
    values_fill = 0
  )

cc_dfs <- split(cc_pivoted, cc_pivoted$Treatment)

cc_matrices <- lapply(cc_dfs, function(df) {
  df %>%
    column_to_rownames("Class") %>%
    .[, -1]                           
})

cc_matrices_clean <- lapply(cc_matrices, function(mat) {
  keep <- apply(mat, 1, function(x) sd(x) > 0 & !all(is.na(x)))
  mat[keep, , drop = FALSE]
})

# keep only taxa present in all groups
common_taxa <- Reduce(intersect, lapply(cc_matrices_clean, rownames))
cc_matrices_clean <- lapply(cc_matrices_clean, function(mat) {
  mat[common_taxa, , drop = FALSE]
})

lapply(cc_matrices_clean, ncol)

cc_k.w_test <- A_star_search(cc_matrices_clean, method = "kw")
cc_perm_test <- A_star_search(cc_matrices_clean, method = "permute")




fert_pivoted <- class_dat %>% 
  filter(Experiment == "Fertilizer") %>% 
  group_by(Class, sample, Treatment) %>%
  summarise(total = sum(Count), .groups = "drop") %>%
  pivot_wider(
    names_from = sample,
    values_from = total,
    values_fill = 0
  )

fert_dfs <- split(fert_pivoted, fert_pivoted$Treatment)

fert_matrices <- lapply(fert_dfs, function(df) {
  df %>%
    column_to_rownames("Class") %>%
    .[, -1]                           
})

zero_var_rows <- unique(unlist(lapply(fert_matrices, function(mat) {
  rownames(mat)[apply(mat, 1, var) == 0]
})))

fert_matrices_clean <- lapply(fert_matrices, function(mat) {
  mat[!(rownames(mat) %in% zero_var_rows), , drop = FALSE]
})


fert_k.w_test <- A_star_search(fert_matrices_clean, method = "kw")
fert_perm_test <- A_star_search(fert_matrices_clean, method = "permute")



HW_LW_modules <- read_rds("capstone modules/HW-LW.rds")

A_star_search(cc_matrices_clean, HW_LW_modules$`16`, method = "kw")



