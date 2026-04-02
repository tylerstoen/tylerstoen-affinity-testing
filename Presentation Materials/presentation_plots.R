library(tidyverse)

# Generate sample data
set.seed(123)

# Group 1: high positive correlation
x1 <- rnorm(100, 0, 1)
y1 <- 0.8 * x1 + rnorm(100, 0, 0.3)

# Group 2: no correlation
x2 <- rnorm(100, 0, 1)
y2 <- rnorm(100, 0, 1)

# Combine into one data frame
df <- data.frame(
  x = c(x1, x2),
  y = c(y1, y2),
  group = rep(c("Condition 1", "Condition 2"), each = 100)
)

ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  facet_wrap(~ group) +
  geom_smooth(method = "lm", se = FALSE, color = "red", size = 1) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Variable X",
    y = "Variable Y"
  )


## 


