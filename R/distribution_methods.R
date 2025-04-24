#' Probability function for IBD Distribution Objects
#'
#' Returns the probability function of a total IBD or segment count distribution.
#'
#' @param x An object of class `total_ibd_dist` or `segment_count_dist`.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric value.
#' @export
d <- function(x, ...){
  UseMethod("d")
}

#' Expectation for IBD Distribution Objects
#'
#' Computes the expectation (mean) of a total IBD or segment count distribution.
#'
#' @param x An object of class `total_ibd_dist` or `segment_count_dist`.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric value.
#' @export
E <- function(x, ...){
  UseMethod("E")
}

#' Variance for IBD Distribution Objects
#'
#' Computes the variance of a total IBD or segment count distribution.
#'
#' @param x An object of class `total_ibd_dist` or `segment_count_dist`.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric value.
#' @export
var <- function(x, ...){
  UseMethod("var")
}

#' @export
var.default <- function(x, ...){
  stats::var(x, ...)
}

#' Standard Deviation for IBD Distribution Objects
#'
#' Computes the standard deviation of a total IBD or segment count distribution.
#'
#' @param x An object of class `total_ibd_dist` or `segment_count_dist`.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric value.
#' @export
#' @export
sd <- function(x, ...){
  UseMethod("sd")
}

#' @export
sd.default <- function(x, ...){
  stats::sd(x, ...)
}

