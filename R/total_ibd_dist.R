#' Compute probability distribution of total IBD
#'
#' The `total_ibd_dist` function computes the probability distribution
#' of the total IBD length (or fraction) over one or more autosomes.
#'
#' For many pedigree relationships, the probability distribution of the
#' total IBD length over one chromosome is a mixture of two point masses
#' (0 and chromosome length) and a continuous density.
#'
#' If `convolve=TRUE` (the default) and `chromosome_length` has length
#' greater than one, the convolution of the distributions will be obtained
#' by FFT using the [`convolve_ibd_dists`] function. Convolution will
#' typically produce a rapidly increasing number of point masses
#' with very small probabilities which are discarded if the
#' probability falls below a threshold of `1e-9`; see [`convolve_ibd_dists`]
#' for details and finer control.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Defaults to [`pedtools::leaves`](pedigree).
#' @param fraction If TRUE, the distribution of the IBD fraction instead of length will be returned. Default is FALSE.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#' @param convolve Should the distribution of the sum (across chromosomes) be obtained?
#' @param ... Additional parameters passed to [`convolve_ibd_dists`] when `convolve=TRUE`.
#' @return object of class `ibd_dist`
#' @examples
#' ## Total IBD and fraction of IBD for a cousin relationship
#' ped_fc <- pedtools::cousinPed()
#'
#' # total IBD length for 100 cM
#' dist_length <- total_ibd_dist(ped_fc, chromosome_length = 100)
#' dist_length
#' plot(dist_length)
#'
#' # fraction IBD for 100 cM
#' dist_fraction <- total_ibd_dist(ped_fc, chromosome_length = 100,
#'                                 fraction = TRUE)
#' dist_fraction
#' plot(dist_fraction)
#'
#' # distribution of total length across three chromosomes (150, 200, 250 cM)
#' plot(total_ibd_dist(ped_fc, chromosome_length = c(150, 200, 250)))
#'
#' # a quick approximation with reasonable accuracy (with just 256 gridpoints)
#' plot(total_ibd_dist(ped_fc, chromosome_length = c(150, 200, 250),
#'                     number_of_gridpoints_exponent = 8))
#'
#' ## Difference between IBD distributions between half-sibs, uncle-nephew
#' ## and grandparent-grandchild relationships
#'
#' # kappa1 is 1/2 for half sibs, uncle-nephew and grandparent-grandchild
#' # but the distribution of the fraction of a chromosome that is in this
#' # state differs between the relationships
#'
#' # define pedigrees and verify kappa1
#' ped_gp <- pedtools::linearPed(n = 2)
#' ped_av <- pedtools::avuncularPed()
#' ped_hs <- pedtools::halfSibPed()
#'
#' stopifnot(all.equal(1/2,
#'           d_ibd(1, ped_av),
#'           d_ibd(1, ped_hs),
#'           d_ibd(1, ped_gp, ids = c(1,5))))
#'
#' # Compute IBD distributions
#' d_av <- total_ibd_dist(ped_av)
#' d_hs <- total_ibd_dist(ped_hs)
#' d_gp <- total_ibd_dist(ped_gp, ids = c(1,5))
#'
#' # the point masses are different
#' d_av
#' d_hs
#' d_gp
#'
#' # plot the continuous densities
#' x0 <- seq(0, 267.77, length.out = 200)
#' df <- data.frame(cM = rep(x0, 3),
#'                  y = c(d(d_av)(x0), d(d_hs)(x0), d(d_gp)(x0)),
#'                  Relationship = rep(c("Avuncular", "Half-Sibling",
#'                                       "Grandparent"), each = length(x0)))
#'
#' require(ggplot2)
#' ggplot(df, aes(x = cM, y = y, color = Relationship)) + geom_line()
#' @importFrom stats setNames qpois dpois
#' @export
total_ibd_dist <- function(pedigree,
                     ids = pedtools::leaves(pedigree),
                     fraction = FALSE,
                     states = "ibd",
                     ibd_state = 1L,
                     chromosome_length = 267.77,
                     convolve = TRUE,
                     ...){

  # validate inputs
  .validate_chromosome_length(chromosome_length)

  states_idx <- .validate_states(states)
  .check_ids_compatible_with_states_idx(ids, states_idx)
  .validate_obs_compatible_with_states_idx(ibd_state, "ibd_state", states_idx)
  .validate_pedigree(pedigree, continuous_genome = TRUE)

  i <- inheritance_space(pedigree = pedigree, ids = ids,
                         states = states)

  lambda_max <- 0.01 * max(chromosome_length) * i$number_of_relevant_transmissions
  joint_n_max <- stats::qpois(1.0 - 1e-16, lambda = lambda_max)

  if (is.infinite(joint_n_max)){
    stop("n_max is infinite")
  }

  # compute the probability distribution of spending k = 0, 1, 2, ..., n_max + 1
  # intervals in the ibd_state if there are n = 0, 1, 2, ..., n_max
  # recombinations
  V <- pr_number_of_intervals_in_state_by_n(ibd_state = ibd_state,
                                            ibd_state_by_v = i$ibd_state_by_v, n_max = joint_n_max,
                                            number_of_transmissions = i$number_of_relevant_transmissions,
                                            masks = i$relevant_masks)
  V_lbeta <- precompute_V_lbeta(V)

  colnames(V) <- paste0("n=", 0:(ncol(V) - 1))
  rownames(V) <- paste0("k=", 0:(nrow(V) - 1))

  fs <- lapply(seq_along(chromosome_length), function(i_chromosome){
    lambda <- 0.01 * chromosome_length[i_chromosome] *
      i$number_of_relevant_transmissions

    n_max <- stats::qpois(1.0 - 1e-16, lambda = lambda)
    n <- 0:n_max
    pr_n <- stats::setNames(stats::dpois(x = n, lambda = lambda), n)

    point_mass_0 <- pr_never_in_state(V = V, n_pr = pr_n, n_max = n_max)
    point_mass_1 <- pr_always_in_state(V = V, n_pr = pr_n, n_max = n_max)

    f_fraction <- function(x){
      d_fraction_ibd_state(s = x, L = chromosome_length[i_chromosome], n_pr = pr_n,
                            V = V, V_lbeta = V_lbeta, point_mass_0, point_mass_1)
    }

    scale <- if (fraction) 1.0 else chromosome_length[i_chromosome]
    f_continuous <- if (fraction) f_fraction else
            function(x) (1/scale) *  f_fraction(x / scale)

    point_mass <- data.frame(x = c(0., scale),
                             px = c(point_mass_0, point_mass_1))
    weight_continuous <- 1 - point_mass_0 - point_mass_1

    dist_chromosome <- list(f_continuous = f_continuous,
                            fraction = fraction,
                            states = states,
                            ibd_state = ibd_state,
                            low = 0,
                            up = scale,
                            chromosome_length = chromosome_length[i_chromosome],
                            weight_continuous = weight_continuous,
                            point_mass = point_mass)

    class(dist_chromosome) <- "ibd_dist"

    dist_chromosome
  })

  if (length(fs) > 1 && convolve){

    number_of_gridpoints_exponent <- 12L
    point_mass_eps <- 1e-9

    args_list <- list(...)
    if (!is.null(args_list$number_of_gridpoints_exponent)){
      number_of_gridpoints_exponent <- args_list$number_of_gridpoints_exponent
    }
    if (!is.null(args_list$point_mass_eps)){
      point_mass_eps <- args_list$point_mass_eps
    }

    return(convolve_ibd_dists(fs, point_mass_eps = point_mass_eps,
                              number_of_gridpoints_exponent = number_of_gridpoints_exponent))
  }

  if (length(chromosome_length) == 1){
    fs <- fs[[1]]
  }

  fs
}

#' @export
d <- function(x, ...){
  UseMethod("d")
}

#' @method d ibd_dist
#' @export
d.ibd_dist <- function(x, ...){
  x$f_continuous
}

#' @export
E <- function(x, ...){
  UseMethod("E")
}

#' @method E ibd_dist
#' @export
#' @importFrom stats integrate
E.ibd_dist <- function(x, m = 1, ...){

  E_mixed <- 0.

  w <- x$weight_continuous
  if (w > 0){
    f <- d(x)

    E_continuous <- stats::integrate(function(x) (x^m) * f(x),
                              lower = x$low, upper = x$up)$val
    E_mixed <- E_mixed + E_continuous * w
  }

  if ((!is.null(x$point_mass)) && (nrow(x$point_mass) > 0)){
    E_discrete <- sum(x$point_mass$x^m * x$point_mass$px)
    E_mixed <- E_mixed + E_discrete
  }

  E_mixed
}

#' @export
var.ibd_dist <- function(x, ...){
  E.ibd_dist(x, m = 2) - E.ibd_dist(x)^2
}

#' @export
sd.ibd_dist <- function(x, ...){
  sqrt(var.ibd_dist(x))
}

#' @importFrom graphics par plot curve grid points axis mtext
#' @export
plot.ibd_dist <- function(x, ...){

  args_list <- list(...)

  if (is.null(args_list$xlim)){
    xlim <- c(min(x$low, x$point_mass$x),
              max(x$up,  x$point_mass$x))
  }
  if (is.null(args_list$n)){
    n <- 200
  }

  w <- x$weight_continuous
  has_continuous_part <- w > 0
  has_discrete_part <- (!is.null(x$point_mass)) && (nrow(x$point_mass) > 0)

  if ((!has_continuous_part) && (!has_discrete_part)){
    stop("distribution has neither a continuous nor a discrete part")
  }

  original_mar <- graphics::par("mar")
  if (has_continuous_part && has_discrete_part){

    # make space for the second y-axis
    graphics::par(mar = c(5, 4, 4, 5))
  }


  if (has_continuous_part){
    f <- d(x)
    graphics::curve(f, from = x$low, to = x$up, n = n, xlim = xlim, ...)
    graphics::grid()
  }

  if (has_continuous_part && has_discrete_part){
    # add second axis
    graphics::par(new = TRUE)

    # add point massses
    graphics::plot(x$point_mass$x, x$point_mass$px,
        xlim = xlim,
         type = "h", lty=2,
         axes = FALSE, xlab="", ylab="")
    graphics::points(x$point_mass$x, x$point_mass$px)

    graphics::axis(4)
    graphics::mtext("p(x)", side = 4, line = 3)

  }else if (has_discrete_part){

    graphics::plot(x$point_mass$x, x$point_mass$px,
         type = "h", ylab="p(x)", lty=2, ...)

    graphics::points(x$point_mass$x, x$point_mass$px)
  }

  graphics::par(mar = original_mar)
}

#' @export
print.ibd_dist = function(x, ...) {
  cat("Probability distribution of",
      if (x$fraction) "fraction" else "total length",
      "of segments in",
      x$states,
      "state",
      x$ibd_state, " \n")
  cat("Chromosome length:", x$chromosome_length, "cM\n\n")
  cat("Weight of continuous density:", x$weight_continuous, "\n\n")

  cat("Point masses: \n")
  print.data.frame(x$point_mass, row.names = FALSE)
}
