#' Probability density of continuously observed IBD segment lengths for pedigree members
#'
#' The `cibd_density` function computes the probability density of observed
#' IBD segments on a chromosome.
#'
#' @param cM Numeric vector with lengths of segments (centiMorgan).
#' @param ibd Integer vector. Taking values 0, 1, 2 for `coefficients = "kappa"`, 1, ..., 9 for `coefficients="identity"` and 1, ..., 15 for `coefficients = "detailed"`.
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param persons Persons for which IBD is observed. Defaults to [`pedtools::leaves`](pedigree).
#' @param coefficients One of `"kappa"`, `"identity"` or `"detailed"`.
#' @param log10 Should the log10 probability density be returned? Default is `FALSE`.
#' @return Numeric
#' @examples
#' ped_fs <- pedtools::nuclearPed(nch = 2)
#'
#' # probability that full siblings are double IBD (kappa2)
#' cibd_density(cM = 0, ibd = 2, ped_fs)
#'
#' # full siblings are double IBD and remain so for 100 cM
#' cibd_density(cM = 100, ibd = 2, ped_fs)
#'
#' # full siblings are double IBD for 50 cM,
#' # then single IBD for 50 cM
#' cibd_density(cM = c(50, 50), ibd = c(2, 1), ped_fs)
#'
#' # full siblings are double IBD, remain so for 100 cM
#' # and no longer
#' cibd_density(cM = c(100, 0), ibd = c(2, 1), ped_fs)
#'
#' ## probability density of IBD segment length for first cousins on an infinite chromosome
#' ped_fc <- pedtools::cousinPed()
#' # first compute the probability of IBD
#' k1_fc <- ibd_pr(ibd = 1, ped_fc)
#' # density of segment length
#' f <- Vectorize(function(x) cibd_density(cM = c(x,0), ibd = c(1, 0), ped_fc) / k1_fc)
#'
#' curve(f, 0, 300)
#'
#' # f is a probability density (integrates to 1)
#' integrate(f, 0, Inf)
#'
#' # for full siblings, how does the chance of remaining double IBD
#' # depend on the segment length?
#' cM <- seq(from = 0, to = 100, length = 200)
#' pr_2ibd <- sapply(cM, cibd_density, 2, ped_fs) / ibd_pr(2, ped_fs)
#'
#' plot(cM, pr_2ibd, type="l")
#'
#' # since there are four meioses, the sojourn time in this IBD state
#' # follows an exponential distribution with rate 0.04
#' stopifnot(all.equal(pr_2ibd, pexp(cM, rate = 0.04, lower.tail = FALSE)))
#' @export
cibd_density <- function(cM, ibd,
                         pedigree, persons = pedtools::leaves(pedigree),
                         coefficients = "kappa",
                         log10 = FALSE){

  # validate inputs
  coeff <- .validate_coefficients(coefficients)
  .check_persons_compatible_with_coeff(persons, coeff)
  .validate_obs_compatible_with_coeff(ibd, "ibd", coeff)
  .validate_logical(log10, "log10")

  i <- inheritance_space(pedigree = pedigree, persons = persons,
                         coefficients = coefficients)

  log10_pr <- log10_ibd_segment_pr_cpp(cM, ibd,
                                       i$ibd_state_by_v,
                                       nrow(i$transmissions),
                                       i$fixed_transmission_masks)

  if (log10){
    return(log10_pr)
  }
  else{
    return(10 ^ log10_pr)
  }
}
