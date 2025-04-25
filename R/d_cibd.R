#' Probability density of continuously observed IBD segment lengths for pedigree members
#'
#' The `d_cibd` function computes the probability density of observed
#' IBD segments on a chromosome.
#'
#' @param x Numeric vector with lengths of segments (centiMorgan) or result from [`r_cibd`].
#' @param ibd Integer vector with IBD states in segments if `x` is not a result from `r_cibd`. Taking values 0, 1, 2 for `states = "ibd"` or `states = "kappa"`, 1, ..., 9 for `states="identity"` and 1, ..., 15 for `states = "detailed"`.
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Default is `pedtools::leaves(pedigree)`.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param log10 Should the log10 probability density be returned? Default is `FALSE`.
#' @return Numeric
#' @examples
#' ped_fs <- pedtools::nuclearPed(nch = 2)
#'
#' # probability that full siblings are double IBD (kappa2)
#' d_cibd(x = 0, ibd = 2, ped_fs)
#'
#' # full siblings are double IBD and remain so for 100 cM
#' d_cibd(x = 100, ibd = 2, ped_fs)
#'
#' # full siblings are double IBD for 50 cM,
#' # then single IBD for 50 cM
#' d_cibd(x = c(50, 50), ibd = c(2, 1), ped_fs)
#'
#' # full siblings are double IBD, remain so for 100 cM
#' # and no longer
#' d_cibd(x = c(100, 0), ibd = c(2, 1), ped_fs)
#'
#' ## probability density of IBD segment length for first cousins on an infinite chromosome
#' ped_fc <- pedtools::cousinPed()
#' # first compute the probability of IBD
#' k1_fc <- d_ibd(ibd = 1, ped_fc)
#' # density of segment length
#' f <- Vectorize(function(l) d_cibd(x = c(l,0), ibd = c(1, 0), ped_fc) / k1_fc)
#'
#' curve(f, 0, 300)
#'
#' # f is a probability density (integrates to 1)
#' integrate(f, 0, Inf)
#'
#' # for full siblings, how does the chance of remaining double IBD
#' # depend on the segment length?
#' cM <- seq(from = 0, to = 100, length = 200)
#' pr_2ibd <- sapply(cM, d_cibd, 2, ped_fs) / d_ibd(2, ped_fs)
#'
#' plot(cM, pr_2ibd, type="l")
#'
#' # since there are four meioses, the sojourn time in this IBD state
#' # follows an exponential distribution with rate 0.04
#' stopifnot(all.equal(pr_2ibd, pexp(cM, rate = 0.04, lower.tail = FALSE)))
#'
#' ## Using the output from r_cibd directly to compute evidential value
#' ## of IBD segments on chromosome for distinguishing first and second cousins
#' ped_c1 <- pedtools::cousinPed()
#' ped_c2 <- pedtools::cousinPed(degree = 2)
#'
#' r_c1 <- r_cibd(n = 1e2, pedigree = ped_c1)
#'
#' lr <- d_cibd(r_c1, pedigree = ped_c1) / d_cibd(r_c1, pedigree = ped_c2)
#'
#' hist(log10(lr))
#
#' @export
d_cibd <- function(x,
                   ibd,
                   pedigree, ids = pedtools::leaves(pedigree),
                   states = "ibd",
                   log10 = FALSE){

  use_rcibd_result <- is.list(x)

  if (use_rcibd_result){
    .validate_rcibd_result(x)

    if (!missing(ibd)){
      stop("when x is a result from r_ibd, the ibd parameter should be missing")
    }

    cM <- x$samples$length
    ibd <- x$samples$state
  }else{
    cM <- x
  }

  # validate inputs
  states_idx <- .validate_states(states)
  .check_ids_compatible_with_states_idx(ids, states_idx)
  .validate_obs_compatible_with_states_idx(ibd, "ibd", states_idx)
  .validate_pedigree(pedigree, continuous_genome = TRUE)
  .validate_logical(log10, "log10")

  i <- inheritance_space(pedigree = pedigree, ids = ids,
                         states = states)

  if (!use_rcibd_result){
    log10_pr <- log10_ibd_segment_pr_cpp(cM, ibd,
                                         i$ibd_state_by_v,
                                         i$number_of_relevant_transmissions,
                                         i$relevant_masks)
  }
  else{
    log10_pr <- log10_ibd_segment_pr_vectorised_cpp(sample = x$samples$sample,
                                        chromosome = x$samples$chromosome,
                                        obs_cM = x$samples$length,
                                        obs_ibd = x$samples$state,
                                        ibd_state_by_v = i$ibd_state_by_v,
                                        number_of_transmissions = i$number_of_relevant_transmissions,
                                        fixed_transmission_masks = i$relevant_masks)
  }

  if (log10){
    return(log10_pr)
  }
  else{
    return(10 ^ log10_pr)
  }
}



