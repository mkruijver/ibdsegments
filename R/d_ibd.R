#' Compute probability of IBD for pedigree members
#'
#' The `d_ibd` function computes the likelihood of IBD
#' for one position or multiple linked markers on the same chromosome.
#'
#' @param ibd Integer vector. Taking values 0, 1, 2 for `states = "ibd"` or `states = "kappa"`, 1, ..., 9 for `states="identity"` and 1, ..., 15 for `states = "detailed"`.
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Defaults to [`pedtools::leaves`]`(pedigree)`.
#' @param recombination_rate_by_locus Numeric vector with length one shorter than `ibd`.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param log10 Should the log10 likelihood be returned? Default is `FALSE`.
#' @return Numeric
#' @examples
#' # Compute kappa0, kappa1, kappa2 for full siblings
#' ped_fs <- pedtools::nuclearPed(nch = 2)
#'
#' k0 <- d_ibd(ibd = 0, pedigree = ped_fs)
#' k1 <- d_ibd(ibd = 1, pedigree = ped_fs)
#' k2 <- d_ibd(ibd = 2, pedigree = ped_fs)
#' c(k0, k1, k2)
#'
#' stopifnot(identical(c(k0, k1, k2), c(0.25, 0.5, 0.25)))
#'
#' # Compute kappa00 for two tightly linked loci
#' d_ibd(c(0,0), pedigree = ped_fs,
#'        recombination_rate_by_locus = c(0.01))
#'
#' # or 100 tightly linked loci
#' d_ibd(rep(0, 100), pedigree = ped_fs,
#'                recombination_rate_by_locus = c(rep(0.01, 99)))
#'
#' # Jacquard's 9 condensed and 15 detailed identity coefficients
#' ped_fs_mating <- pedtools::fullSibMating(1)
#'
#' sapply(1:9, d_ibd, pedigree = ped_fs_mating, states = "identity")
#' sapply(1:15, d_ibd, pedigree = ped_fs_mating, states = "detailed")
#' @export
d_ibd <- function(ibd,
                   pedigree, ids = pedtools::leaves(pedigree),
                   recombination_rate_by_locus = numeric(),
                   states = "ibd",
                   log10 = FALSE){

  states_idx <- .validate_states(states)
  .check_ids_compatible_with_states_idx(ids, states_idx)
  .validate_obs_compatible_with_states_idx(ibd, "ibd", states_idx)
  validate_recombination_rates_cpp(recombination_rate_by_locus)
  .validate_recombination_rates_compatible_with_obs(ibd,
    "recombination_rate_by_locus", recombination_rate_by_locus)
  .validate_logical(log10, "log10")

  if (length(recombination_rate_by_locus) > 0){
    .assert_no_founder_inbreeding(pedigree,
        "Founder inbreeding is not supported for calculations involving more than one locus")
  }


  i <- inheritance_space(pedigree = pedigree, ids = ids,
                         states = states)

  if (pedtools::hasInbredFounders(pedigree)){
    .assert_ids_are_not_inbred_founders(pedigree, ids)

    pr_v_constant <- -1
    pr_v <- v_prior_with_f(pedigree, i)
  }else{
    pr_v_constant <- 2^(-(i$number_of_relevant_transmissions - length(i$relevant_masks)))
    pr_v <- numeric()
  }

  log10_likelihood <- sum(ibd_log10_pr_cpp(ibd_state_by_v = i$ibd_state_by_v,
                                            ibd_by_locus = ibd,
                                            recombination_rate_by_locus = recombination_rate_by_locus,
                                            number_of_transmissions = i$number_of_relevant_transmissions,
                                            fixed_transmission_masks = i$relevant_masks,
                                           pr_v_constant = pr_v_constant,
                                           pr_v_biased = pr_v))

  if (log10){
    return(log10_likelihood)
  }
  else{
    return(10 ^ log10_likelihood)
  }
}
