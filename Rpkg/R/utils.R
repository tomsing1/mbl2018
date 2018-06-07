#' Tests if a vector is "categorical"
#'
#' Is this thing a character or factor?
#'
#' @export
#' @param x a vector to test
#' @return \code{TRUE} if \code{x} is categorical
is.categorical <- function(x) {
  is.character(x) || is.factor(x)
}

#' Tests if a vector is "integerish"
#'
#' Is this thing integer like? As in, is it an integer, or a numeric that has
#' no decimals?
#'
#' @export
#' @param x a vector to test
#' @return \code{TRUE} if \code{x} is integerish.
is.integerish <- function(x) {
  test_integerish(x)
}
