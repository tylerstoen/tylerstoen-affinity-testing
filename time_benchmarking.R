library(bench)
library(tidyverse)

# Time Checks for number of rows.
results <- bench::press(
  n = seq(10, 100, 10),
  {
    bench::mark(
      {
        test_data <- test_maker(n, 20, 0.7, 0.2, 0.1, c(1:10))
        find_A_star_delta2(test_data, c("V1", "V2", "V5"))
      },
      iterations = 10,
      check = FALSE
    )
  }
)

results

results |>
  mutate(
    median_time = as.double(median)
  ) |>
  ggplot() +
  aes(x = n, y = median_time) +
  geom_point() +
  geom_line() +
  theme_bw()



