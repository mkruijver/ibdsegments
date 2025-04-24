#' Convolve IBD distributions to obtain the distribution of the sum
#'
#' Chromosome-wise IBD distributions obtained from the `total_ibd_dist`
#' function can be convoluted manually using the `convolve_ibd_dists`
#' function. This allows finer control of the procedure by controlling
#' the number of gridpoints used in the FFT and the threshold for
#' point masses to be retained.
#'
#' The `convolve_ibd_dist` implements a convolution procedure based on a self
#' contained variant of the recipe implemented in the `distr` package adapted
#' to this use case. In particular, the IBD distribution for one
#' chromosome is often a mixture of two point masses and continuous with finite
#' support otherwise. Convolution increases the number of point masses and
#' decreases their probability mass. This function assumes that small point
#' masses are not of specific interest and will discard the point masses when
#' the probability mass is smaller than `point_mass_eps` with a default value
#' of `1e-9`. For typical pedigree relationships this means that the IBD
#' distribution of more than a few chromosomes is treated as a continuous
#' distribution admitting a density function.
#'
#' The `number_of_gridpoints_exponent` controls the density of gridpoints in the
#' FFT. By default this is `12` which means that `2^12=4096` gridpoints are
#' used. Increasing this parameter improves accuracy at the cost of increased
#' runtime and memory use.
#'
#' The point masses are of particular interest in some applications. For
#' example, the probability that no part of the autosomal genome is inherited
#' may be of interest. In that case, the point masses should not be discarded
#' and the `d_cibd` function may be used. See the example below, for further
#' details.
#'
#' @param ... ibd dists.
#' @param number_of_gridpoints_exponent Default is 12.
#' @return ibd_dist
#' @examples
#' ## Convolution of IBD distributions for half siblings at chromosome 1 and 2
#'
#' # define half sib pedigree
#' ped_hs <- pedtools::halfSibPed()
#'
#' # obtain chromosome-wise distributions
#' d1 <- total_ibd_dist(pedigree = ped_hs, chromosome_length = 267.77)
#' d2 <- total_ibd_dist(pedigree = ped_hs, chromosome_length = 251.73)
#'
#' convolve_ibd_dists(d1, d2) # 4 point masses
#'
#' ## Accuracy of convolution depends on number of gridpoints
#'
#' L <- c(267.77, 251.73, 218.31, 202.89, 197.08, 186.02, 178.4, 161.54,
#'        157.35, 169.28, 154.5, 165.49, 127.23, 116, 117.32, 126.59, 129.53,
#'        116.52, 106.35, 107.76, 62.88, 70.84)
#'
#' # obtain chromosome-wise distributions for half siblings
#' hs <- total_ibd_dist(pedigree = ped_hs,
#'                      chromosome_length = L, convolve = FALSE)
#'
#' # obtain moments of the sum by summing the moments..
#' mean_hat <- sum(sapply(hs, E))
#' sd_hat <- sqrt(sum(sapply(hs, var)))
#'
#' # .. or by summing the distributions with varying numbers of gridpoints
#' k <- 6:10
#' sd_hat_by_k <- sapply(k, function(k) sd(convolve_ibd_dists(hs,
#'                                         number_of_gridpoints_exponent = k)))
#' mean_hat_by_k <- sapply(k, function(k) E(convolve_ibd_dists(hs,
#'                                         number_of_gridpoints_exponent = k)))
#'
#' plot(k, mean_hat_by_k)
#' abline(h = mean_hat, lty = 2)
#'
#' plot(k, sd_hat_by_k)
#' abline(h = sd_hat, lty = 2)
#' @export
convolve_ibd_dists <- function(...,
                               point_mass_eps = 1e-9,
                               number_of_gridpoints_exponent = 12){
  x <- list(...)
  if (length(x) == 1){
    if (is.list(x[[1]])){
      x <- x[[1]]
    }
  }

  # start with the smallest and work your way up
  scales <- sapply(x, function(d0) d0$up)
  o <- order(scales)

  partial_sum <- x[[o[1]]]

  for(i in o[-1]){
    partial_sum <- convolve_two_ibd_dists(partial_sum,
                                          x[[i]], point_mass_eps,
                                          number_of_gridpoints_exponent)
  }

  partial_sum
}

#' @importFrom stats fft approxfun
convolve_two_ibd_dists <- function(d1, d2, point_mass_eps, number_of_gridpoints_exponent){
  .validate_ibs_dist(d1)
  .validate_ibs_dist(d2)
  .validate_not_ibd_fraction_dist(d1)
  .validate_not_ibd_fraction_dist(d2)

  n_gridpoints <- 2^number_of_gridpoints_exponent

  # step 1: truncation
  new_lower <- min(d1$low, d2$low)
  new_upper <- max(d1$up, d2$up)
  h <- (new_upper - new_lower) / n_gridpoints

  # step 2: discretise both distributions on the same grid
  f1 <- discretise_pdf(f = d1$f_continuous, low = d1$low, up = d1$up,
                       n_gridpoints = n_gridpoints)
  f2 <- discretise_pdf(f = d2$f_continuous, low = d2$low, up = d2$up,
                       n_gridpoints = n_gridpoints)

  p1 <- integrate_pdf_to_cdf(f = f1, low = d1$low,
                             up = d1$up, n_gridpoints = n_gridpoints)

  p2 <- integrate_pdf_to_cdf(f = f2, low = d2$low,
                             up = d2$up, n_gridpoints = n_gridpoints)

  dp1 <- discretise_cdf(p = p1, p_lower = d1$low, p_upper = d1$up,
                        new_lower = new_lower, new_upper = new_upper,
                        h = h, n_gridpoints)
  dp2 <- discretise_cdf(p = p2, p_lower = d2$low, p_upper = d2$up,
                        new_lower = new_lower, new_upper = new_upper,
                        h = h, n_gridpoints)

  # step 3: zero padding
  dp1_padded <- c(dp1, numeric(n_gridpoints))
  dp2_padded <- c(dp2, numeric(n_gridpoints))

  # step 4: transform
  f1_hat <- stats::fft(dp1_padded)
  f2_hat <- stats::fft(dp2_padded)

  # convolution in Fourier domain
  f12_hat <- f1_hat * f2_hat

  d_sum_cont <- c(0, Re(stats::fft(f12_hat, inverse = TRUE)) / length(f1_hat)) / h

  # blitz out artefacts
  d_sum_cont[d_sum_cont < .Machine$double.eps] <- 0.
  x_grid <- seq(from = 2*new_lower,
                to = 2 * new_upper, by = h)

  point_mass <- pmf_of_sum(x1 = d1$point_mass$x, p1 = d1$point_mass$px,
                           x2 = d2$point_mass$x, p2 = d2$point_mass$px,
                           eps = point_mass_eps)

  sum_weight_cont <- 1 - sum(point_mass$px)

  sum_weight_cc <- d1$weight_cont * d2$weight_cont
  sum_weight_pc <- sum(d1$point_mass$px) * d2$weight_cont
  sum_weight_cp <- sum(d2$point_mass$px) * d1$weight_cont

  # the continuous part of the distribution of the sum
  # comprises three components:
  # both summands continuous, discrete + cont and cont + discrete
  sum_weights <- sum_weight_cc  + sum_weight_pc + sum_weight_cp

  # sum_weight_pc
  d_sum <- d_sum_cont * sum_weight_cc/sum_weights
  for (i_point in seq_len(nrow(d1$point_mass))){
    x_point <- d1$point_mass$x[i_point]
    px_point <- d1$point_mass$px[i_point]

    d_sum <- d_sum + (px_point * d2$weight_cont)/sum_weights * f2(x_grid - x_point)
  }

  for (i_point in seq_len(nrow(d2$point_mass))){
    x_point <- d2$point_mass$x[i_point]
    px_point <- d2$point_mass$px[i_point]

    d_sum <- d_sum + (px_point * d1$weight_cont)/sum_weights * f1(x_grid - x_point)
  }

  sum_lower <- x_grid[find_index_of_first_non_zero(d_sum, eps = .Machine$double.eps)]
  sum_upper <- x_grid[find_index_of_last_non_zero(d_sum, eps = .Machine$double.eps)]

  f_sum <- stats::approxfun(x = x_grid, y = d_sum, yleft = 0, yright = 0)

  new_length <-  c(d1$chromosome_length,
                   d2$chromosome_length)

  dist_sum <- list(f_continuous = f_sum,
                   fraction = d1$fraction,
                   states = d1$states,
                   ibd_state = d1$ibd_state,
                   low = sum_lower,
                   up = sum_upper,
                   chromosome_length = new_length,
                   weight_continuous = sum_weight_cont,
                   point_mass = point_mass)
  class(dist_sum) <- "ibd_dist"

  dist_sum
}

discretise_pdf <- function(f, low, up, n_gridpoints){
  xx <- seq(low, up, length = n_gridpoints + 1)
  dx <- f(xx)

  approxfun(x = xx, y = dx, yleft = 0, yright = 0)
}

#' @importFrom stats approxfun
integrate_pdf_to_cdf <- function(f, low, up, n_gridpoints){
  # this is distr:::.D2P

  xx <- seq(low, up, length = n_gridpoints + 1)
  dx <- f(xx)

  # h_p <- (up - low) / (n_gridpoints + 1)

  x_less_than <- xx[seq_along(xx) %% 2 == 1]
  p_less_than <- cumulative_simpson_cpp(dx)
  p_less_than <- p_less_than / max(p_less_than)

  stats::approxfun(x = x_less_than, y = p_less_than,
            yleft = 0, yright = 1)
}

discretise_cdf <- function(p, p_lower, p_upper, new_lower, new_upper, h, n_gridpoints){
  # this is distr:::..discretizeP

  h0 <- 40 * (p_upper - p_lower) * (1 / n_gridpoints)

  if (h > h0)
    warning("Grid for approxfun too wide, increase number of grid points")

  x <- seq(from = new_lower, to = new_upper, by = h)

  return(diff(p(x)))
}
