library(dplyr)
example_dge_data <- function(seed = 300) {
  set.seed(seed)
  tibble(
    gene = rep(c("Gene A", "Gene B", "Gene C"), each = 6),
    condition = factor(rep(rep(c("WT", "KO"), each = 3), times = 3), c("WT" ,"KO")),
    expression = c(
      rnorm(3, mean = 2, sd = 1), rnorm(3, mean = -3, sd = 1),     # Gene A
      rnorm(3, mean = 1, sd = 2), rnorm(3, mean = -1, sd = 2),     # Gene B
      rnorm(3, mean = 0, sd = 1.25), rnorm(3, mean = 0, sd = 1.25) # Gene C
    ))
}

