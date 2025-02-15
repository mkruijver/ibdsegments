#' Compute probability distribution of total IBD
#'
#' The `ibd_dist` function computes the probability distribution
#' of the total IBD fraction over an autosome. A function is returned
#' that gives the probability density at 0 < x < 1 and the point masses
#' at x=0 and x=1.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param persons Persons for which IBD is observed. Defaults to [`pedtools::leaves`](pedigree).
#' @param fraction If TRUE, the distribution of the IBD fraction instead of length will be returned. Default is FALSE.
#' @param coefficients One of `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#' @param convolve Should the distribution of the sum (across chromosomes) be obtained?
#' @return object of class `ibd_dist`
#' @examples
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
#'           ibd_pr(1, ped_av),
#'           ibd_pr(1, ped_hs),
#'           ibd_pr(1, ped_gp, persons = c(1,5))))
#'
#' # Compute IBD distributions
#' d_av <- ibd_dist(ped_av)
#' d_hs <- ibd_dist(ped_hs)
#' d_gp <- ibd_dist(ped_gp, persons = c(1,5))
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
#' @export
ibd_dist <- function(pedigree,
                     persons = pedtools::leaves(pedigree),
                     fraction = FALSE,
                     coefficients = "kappa",
                     ibd_state = 1L,
                     chromosome_length = 267.77,
                     convolve = TRUE){

  # validate inputs
  if (!is.numeric(chromosome_length)){
    stop("Chromosome_length needs to be numeric")
  }

  if (any(chromosome_length <= 0)){
    stop("Chromosome_length needs to be strictly positive")
  }

  coeff <- .validate_coefficients(coefficients)
  .check_persons_compatible_with_coeff(persons, coeff)
  .validate_obs_compatible_with_coeff(ibd_state, "ibd_state", coeff)

  i <- inheritance_space(pedigree = pedigree, persons = persons,
                         coefficients = coefficients)

  lambda_max <- 0.01 * max(chromosome_length) * i$number_of_relevant_transmissions
  joint_n_max <- qpois(1.0 - 1e-16, lambda = lambda_max)

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

    n_max <- qpois(1.0 - 1e-16, lambda = lambda)
    n <- 0:n_max
    pr_n <- setNames(dpois(x = n, lambda = lambda), n)

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
                            coefficients = coefficients,
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
    return(convolve_ibd_dists(fs))
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

#' @export
d.ibd_dist <- function(x){
  x$f_continuous
}

#' @export
print.ibd_dist = function(x, ...) {
  cat("Probability distribution of",
      if (x$fraction) "fraction" else "total length",
      "of segments in",
      x$coefficients,
      "state",
      x$ibd_state, " \n")
  cat("Chromosome length:", x$chromosome_length, "cM\n\n")
  cat("Weight of continuous density:", x$weight_continuous, "\n\n")

  cat("Point masses: \n")
  print.data.frame(x$point_mass, row.names = FALSE)
}
