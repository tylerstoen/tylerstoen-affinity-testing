library(here)
library(tidyverse)
library(paletteer)
library(broom)
library(vegan)
library(ggplot2)
library(patchwork)
library(stringr)

get_correlation_matrix <- function(data, taxa, treatment){
  
  cb_trt <- data |>
    filter(Treatment == treatment) |>
    select(all_of(taxa), Sample_ID, Count) |>
    pivot_wider(names_from = Sample_ID, values_from = Count, values_fill = 0)
  
  taxa_labels <- cb_trt[[taxa]]
  
  
  use <- cb_trt |>
    select(-Class)
  
  filtered <- use[, apply(use, 2, sd) != 0]
  
  # Compute the correlation matrix if no constant species remain
  correlation_matrix <- cor(t(filtered))
  
  colnames(correlation_matrix) <- taxa_labels
  rownames(correlation_matrix) <- taxa_labels
  
  return(correlation_matrix)
}


get_correlation_heatmap_2 <- function(matrix, taxa){
  heatmap <- matrix[taxa, taxa]
  heatmap[is.na(heatmap)] <- 0
  library(reshape2)
  
  heatmap_long <- melt(heatmap)
  
  ggplot(heatmap_long, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 0)) +
    xlab("") + 
    ylab("") +
    scale_y_discrete(limits=rev) +
    scale_x_discrete(position = "top")
}
