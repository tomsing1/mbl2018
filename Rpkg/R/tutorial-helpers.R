#' Generates example data for three categories of differential expression
#'
#' @export
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

#' Simulation of data measured for "mouse tail length experiment"
#'
#' @description
#' This function simulates measurements of tail lengths from two different
#' groups of mice. Mice in the `"control"` group were given a normal diet, and
#' mice in the `"treatment"` group were given a diet with a supplement that
#' claims to increase tail length.
#'
#' Explore how the number of mice measured per group effects the conclusion from
#' the experiment by playing with (1) the number of mice used per experimental
#' group; and (2) the *actual* difference tail length in the treatment group.
#'
#' This function assumes that the average tail length of a mouse is 8cm.
#'
#' @export
#' @param mice_per_group The number of mice used per group
#' @param actual_increase The actual increase in length (cm) the supplement
#'   has on the tail length of mice. The default value (1) indicates
#'   **an actuall** increase of 1cm in tail length of the treatment mice.
#' @return a data.frame of tail length measurements for the mice in the
#'   experiment.
run_tail_length_experiment <- function(num_mice = 3, actual_increase = 1,
                                       variance = 1) {
  # variance is square of sd
  ctrl.length <- rnorm(num_mice, mean = 6, sd = sqrt(variance))
  trt.length <- rnorm(num_mice, mean = 6 + actual_increase, sd = sqrt(variance))

  data <- data.frame(
    group = rep(c("control", "treatment"), each = num_mice),
    tail_length = c(ctrl.length, trt.length),
    mouse_id = seq(2*num_mice))
  data
}

if (FALSE) {
  pvals <- sapply(1:1000, function(i) {
    dat <- run_tail_length_experiment(5, 0)
    tt <- t.test(tail_length ~ group, dat)
    tt$p.value
  })

  hist(pvals, 30)
}
