#' Compute probabilities of minimal multi-person IBD patterns
#'
#' For two full siblings at one locus, it is well known that there are three distinct minimal IBD patterns
#' with probabilities 0.25, 0.5 and 0.25. The [`ribd::multiPersonIBD`] function
#' generalises the computation of these patterns and their probabilities to more than
#' two ids The `multi_ibd_patterns` function further generalises the computation
#' to patterns across multiple loci.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Defaults to [`pedtools::leaves`]`(pedigree)`.
#' @param recombination_rate_by_locus Optionally a numeric vector with recombination rates.
#' @return DataFrame
#' @examples
#' # Compute IBD patterns for two full siblings...
#' multi_ibd_patterns(pedtools::nuclearPed(nch = 2))
#'
#' # ... and the generalisation to three siblings
#' multi_ibd_patterns(pedtools::nuclearPed(nch = 3))
#'
#' # Two full siblings at two tightly linked loci
#' multi_ibd_patterns(pedtools::nuclearPed(nch = 2),
#'                    recombination_rate_by_locus = 0.01)
#' @export
multi_ibd_patterns <- function(pedigree, ids = pedtools::leaves(pedigree),
                               recombination_rate_by_locus = numeric()){


  # validate inputs
  person_idx <- .validate_pedigree_ids(ids, pedigree)
  validate_recombination_rates_cpp(recombination_rate_by_locus)

  # pedigree will be validated here
  i <- inheritance_space(pedigree, exploit_symmetries = TRUE)

  # generate all minimal multi person IBD patterns by v in the columns of m
  m <- get_multi_ibd_patterns_by_v(number_of_ped_members = length(pedigree$ID),
                                   ped_row_is_founder_idx = which(pedigree$FIDX == 0),
                                   from_allele_idx = i$transmissions$from_allele_idx,
                                   to_allele_idx = i$transmissions$to_allele_idx,
                                   person_ids = person_idx,
                                   number_of_fixed_transmissions = sum(i$transmissions$is_fixed),
                                   top_to_bottom_order = i$transmissions$top_to_bottom_order,
                                   minimal = TRUE)

  number_of_loci <- 1 + length(recombination_rate_by_locus)

  if (number_of_loci == 1){

    if (pedtools::hasInbredFounders(pedigree)){
      .assert_persons_are_not_inbred_founders(pedigree, ids)

      pr_v_constant <- -1
      pr_v <- v_prior_with_f(pedigree, i)
    }else{
      pr_v_constant <- 2^(-(i$number_of_relevant_transmissions - length(i$relevant_masks)))
      pr_v <- numeric()
    }

    return(multi_ibd_patterns_by_v_df(m$unique_patterns, pattern_idx_by_v = m$pattern_idx_by_v,
                                      ids = ids, pr_v_constant = pr_v_constant, pr_v = pr_v))
  }
  else{
    .assert_no_founder_inbreeding(pedigree,
        "Founder inbreeding is not supported for calculations involving more than one locus")

    number_of_states <- max(m$pattern_idx_by_v)

    multi_locus_m_idxs <- as.matrix(rev(do.call(expand.grid,
      replicate(n = number_of_loci, seq_len(number_of_states), simplify = FALSE))))

    prob <- apply(multi_locus_m_idxs, 1, function(multi_locus_m_idx){
        10^sum(ibd_log10_pr_cpp(ibd_state_by_v = m$pattern_idx_by_v,
                             ibd_by_locus = multi_locus_m_idx,
                             recombination_rate_by_locus = recombination_rate_by_locus,
                             number_of_transmissions = i$number_of_relevant_transmissions,
                             fixed_transmission_masks = i$relevant_masks))})

    return(multi_ibd_patterns_df(prob = prob, multi_locus_m_idx = multi_locus_m_idxs,
                          unique_patterns = m$unique_patterns, ids = ids))
  }
}
