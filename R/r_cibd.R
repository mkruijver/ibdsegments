#' Random generation for IBD on a continuous genome
#'
#' The `r_cibd` function generates random Identity-by-Descent (IBD) segments
#' along a continuous genome, given a specified pedigree and observed
#' individuals.
#'
#' @param n Number of observations
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Default is `pedtools::leaves(pedigree)`.
#' @param states One of `"ibd` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#'
#' @return A list containing:
#'   \item{samples}{Data frame of simulated IBD segments along the chromosome.}
#'   \item{stats}{Data frame with summary statistics per sample,
#'   including total IBD length and the segment count.}
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
#' r_gp <- r_cibd(n = 1e4, pedigree = pedtools::linearPed(2), ids = c(1, 5),
#'                chromosome_length = 300)
#'
#' hist(r_gp$stats$total_length)
#' hist(r_hs$stats$total_length)
#' @export
r_cibd <- function(n,
                   pedigree,
                   ids = pedtools::leaves(pedigree),
                   states = "ibd",
                   ibd_state = 1L,
                   chromosome_length = 267.77){

  states_idx <- .validate_states(states)
  .check_ids_compatible_with_states_idx(ids, states_idx)
  .validate_obs_compatible_with_states_idx(ibd_state, "ibd_state", states_idx)
  .validate_pedigree(pedigree, continuous_genome = TRUE)

  i <- inheritance_space(pedigree = pedigree, ids = ids,
                         states = states)

  random_ibd(n = n, chromosome_length = chromosome_length,
             ibd_state_by_v = i$ibd_state_by_v,
             number_of_transmissions = i$number_of_relevant_transmissions,
             fixed_transmission_masks = i$relevant_masks,
             state_stats = ibd_state)
}
