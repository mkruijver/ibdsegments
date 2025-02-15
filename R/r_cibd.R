#' Random generation for IBD on a continuous genome
#'
#' The `r_cibd`
#'
#' @param n Number of observations
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param persons Persons for which IBD is observed. Defaults to [`pedtools::leaves`](pedigree).
#' @param coefficients One of `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#'
#' @examples
#' ## Basic example: IBD along one chromosome for half siblings
#' L <- 300
#' r_hs <- r_cibd(n = 1e4, pedigree = pedtools::halfSibPed(), chromosome_length = 300)
#'
#' # half sibs alternate between IBD (state 1) and not IBD (state 0)
#' head(r_hs$samples)
#'
#' # the total_length and number of segments per sample are also returned
#' head(r_hs$stats)
#'
#' ## Comparing half siblings and grandparent-grandchild
#' r_gp <- r_cibd(n = 1e4, pedigree = pedtools::linearPed(2), persons = c(1, 5),
#'                chromosome_length = 300)
#'
#' hist(r_gp$stats$total_length)
#' hist(r_hs$stats$total_length)
#' @export
r_cibd <- function(n,
                   pedigree,
                   persons = pedtools::leaves(pedigree),
                   coefficients = "kappa",
                   ibd_state = 1L,
                   chromosome_length = 267.77){

  coeff <- ibdsegments::: .validate_coefficients(coefficients)
  .check_persons_compatible_with_coeff(persons, coeff)
  .validate_obs_compatible_with_coeff(ibd_state, "ibd_state", coeff)

  i <- inheritance_space(pedigree = pedigree, persons = persons,
                         coefficients = coefficients)

  random_ibd(n = n, chromosome_length = chromosome_length,
             ibd_state_by_v = i$ibd_state_by_v,
             number_of_transmissions = i$number_of_relevant_transmissions,
             fixed_transmission_masks = i$relevant_masks,
             state_stats = ibd_state)
}
