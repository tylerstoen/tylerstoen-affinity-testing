library(here)
library(tidyverse)
raw_class_counts <- read.csv("rarified_counts_by_Class.csv")

zeros_CC <- raw_class_counts |>
  filter(Experiment == "Crop Cover") |>
  filter(Count == 0) |>
  select(Class) |>
  unique()

filtered_class_counts <- raw_class_counts |>
  filter(Experiment == "Crop Cover") |>
  filter(!Class %in% zeros_CC$Class)

n_distinct(filtered_class_counts$Class)




CC_corr_matrices <- list(get_correlation_matrix(filtered_class_counts, "Class", "HW"),
                         get_correlation_matrix(filtered_class_counts, "Class", "LW"),
                         get_correlation_matrix(filtered_class_counts, "Class", "NCC"))


test <- compute_delta(CC_corr_matrices, "alphaproteobacteria", c("betaproteobacteria"))

print(test)

HW_corr_matrix <- get_correlation_matrix(filtered_class_counts, "Class", "HW")
LW_corr_matrix <- get_correlation_matrix(filtered_class_counts, "Class", "LW")
NCC_corr_matrix <- get_correlation_matrix(filtered_class_counts, "Class", "NCC")
first5 <- rownames(HW_corr_matrix)[1:5]
second5 <- rownames(HW_corr_matrix)[6:10]

get_correlation_heatmap_2(HW_corr_matrix, second5)
get_correlation_heatmap_2(LW_corr_matrix, second5)
get_correlation_heatmap_2(NCC_corr_matrix, second5)

test2 <- compute_delta(CC_corr_matrices, "deltaproteobacteria", c("acidobacteriia"))

test
test2

test_astar_og <- find_A_star_delta(CC_corr_matrices, c("cloacimonetes", "gemmatimonadetes"))

get_correlation_heatmap_2(HW_corr_matrix, test_astar_og)
get_correlation_heatmap_2(LW_corr_matrix, test_astar_og)
get_correlation_heatmap_2(NCC_corr_matrix, test_astar_og)

