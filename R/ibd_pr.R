#' Compute probability of IBD for pedigree members
#'
#' The `ibd_pr` function computes the likelihood of IBD
#' for one position or multiple linked markers on the same chromosome.
#'
#' @param ibd Integer vector. Taking values 0, 1, 2 for `coefficients = "ibd"` or `coefficients = "kappa"`, 1, ..., 9 for `coefficients="identity"` and 1, ..., 15 for `coefficients = "detailed"`.
#' @param pedigree Pedigree in [`pedtools::ped`] form.kappa
#' @param persons Persons for which IBD is observed. Defaults to [`pedtools::leaves`]`(pedigree)`.
#' @param recombination_rate_by_locus Numeric vector with length one shorter than `ibd`.
#' @param coefficients One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param log10 Should the log10 likelihood be returned? Default is `FALSE`.
#' @return Numeric
#' @examples
#' # Compute kappa0, kappa1, kappa2 for full siblings
#' ped_fs <- pedtools::nuclearPed(nch = 2)
#'
#' k0 <- ibd_pr(ibd = 0, pedigree = ped_fs)
#' k1 <- ibd_pr(ibd = 1, pedigree = ped_fs)
#' k2 <- ibd_pr(ibd = 2, pedigree = ped_fs)
#' c(k0, k1, k2)
#'
#' stopifnot(identical(c(k0, k1, k2), c(0.25, 0.5, 0.25)))
#'
#' # Compute kappa00 for two tightly linked loci
#' ibd_pr(c(0,0), pedigree = ped_fs,
#'        recombination_rate_by_locus = c(0.01))
#'
#' # or 100 tightly linked loci
#' ibd_pr(rep(0, 100), pedigree = ped_fs,
#'                recombination_rate_by_locus = c(rep(0.01, 99)))
#'
#' # Jacquard's 9 condensed and 15 detailed identity coefficients
#' ped_fs_mating <- pedtools::fullSibMating(1)
#'
#' sapply(1:9, ibd_pr, pedigree = ped_fs_mating, coefficients = "identity")
#' sapply(1:15, ibd_pr, pedigree = ped_fs_mating, coefficients = "detailed")
#' @export
ibd_pr <- function(ibd,
                   pedigree, persons = pedtools::leaves(pedigree),
                   recombination_rate_by_locus = numeric(),
                   coefficients = "ibd",
                   log10 = FALSE){

  coeff <- .validate_coefficients(coefficients)
  .check_persons_compatible_with_coeff(persons, coeff)
  .validate_obs_compatible_with_coeff(ibd, "ibd", coeff)
  validate_recombination_rates_cpp(recombination_rate_by_locus)
  .validate_recombination_rates_compatible_with_obs(ibd,
    "recombination_rate_by_locus", recombination_rate_by_locus)
  .validate_logical(log10, "log10")

  i <- inheritance_space(pedigree = pedigree, persons = persons,
                         coefficients = coefficients)

  log10_likelihood <- sum(ibd_log10_pr_cpp(ibd_state_by_v = i$ibd_state_by_v,
                                            ibd_by_locus = ibd,
                                            recombination_rate_by_locus = recombination_rate_by_locus,
                                            number_of_transmissions = i$number_of_relevant_transmissions,
                                            fixed_transmission_masks = i$relevant_masks))

  if (log10){
    return(log10_likelihood)
  }
  else{
    return(10 ^ log10_likelihood)
  }
}
