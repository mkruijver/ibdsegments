#' Compute moments of probability distribution of total IBD
#'
#' The `total_ibd_dist_moments` function computes mean and variance of
#' the probability distribution of the total IBD length (or fraction)
#' over one autosome. The function uses double numerical integration.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Defaults to [`pedtools::leaves`](pedigree).
#' @param fraction If TRUE, the distribution of the IBD fraction instead of length will be returned. Default is FALSE.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#' @return `list`
#' @examples
#' # Full Siblings are double IBD with 25% probability
#' # we may compute the expectation and variance of the total length of
#' # a chromosome that is double IBD for a chromosome of 100 cM (i.e. 1 Morgan)
#' m <- total_ibd_dist_moments(pedigree = pedtools::nuclearPed(nch = 2),
#'                        ibd_state = 2, chromosome_length = 100)
#' m
#'
#' # compare to numerical integration of the full distribution
#' d_fs <- total_ibd_dist(pedigree = pedtools::nuclearPed(nch = 2),
#' ibd_state = 2, chromosome_length = 100)
#'
#' m2 <- list(mean = E(d_fs), variance = var(d_fs), sd = sd(d_fs))
#' m2
#'
#' stopifnot(all.equal(m, m2))
#'
#' # Expectation and variance of _fraction_ of the genome that is
#' # double IBD between four full siblings
#' m4 <- total_ibd_dist_moments(pedigree = pedtools::nuclearPed(nch = 4),
#' ibd_state = 2, chromosome_length = 100, fraction = TRUE)
#' m4
#'
#' stopifnot(all.equal(0.25^3, m4$mean))
#' @export
total_ibd_dist_moments <- function(pedigree,
                                     ids = pedtools::leaves(pedigree),
                                     fraction = FALSE,
                                     states = "ibd",
                                     ibd_state = 1L,
                                     chromosome_length = 267.77){

  # validate inputs
  if (!is.numeric(chromosome_length)){
    stop("Chromosome_length needs to be numeric")
  }

  if (any(chromosome_length <= 0)){
    stop("Chromosome_length needs to be strictly positive")
  }
  if (length(chromosome_length) != 1){
    stop("Chromosome_length needs to have length 1")
  }


  states_idx <- .validate_states(states)
  .check_ids_compatible_with_states_idx(ids, states_idx)
  .validate_obs_compatible_with_states_idx(ibd_state, "ibd_state", states_idx)
  .validate_pedigree(pedigree, continuous_genome = TRUE)

  i <- inheritance_space(pedigree = pedigree, ids = ids,
                         states = states)


  # expectation of square of indicator function
  pr_v_constant <- 2^(-(i$number_of_relevant_transmissions - length(i$relevant_masks)))
  pr_v <- numeric()

  ii <- rep(ibd_state, 2)

  E_I_squared <- 1 / (0.5 * chromosome_length^2) *
    integrate(f = Vectorize(function(x){
      integrate(f = Vectorize(function(y) {

        rho <-  0.5 * (1 - exp(-2 * 0.01* (x - y)))

        10^sum(ibd_log10_pr_cpp(ibd_state_by_v = i$ibd_state_by_v,
                             ibd_by_locus = ii,
                             recombination_rate_by_locus = rho,
                             number_of_transmissions = i$number_of_relevant_transmissions,
                             fixed_transmission_masks = i$relevant_masks,
                             pr_v_constant = pr_v_constant,
                             pr_v_biased = pr_v))
      }), lower = 0, upper = x)$val
    }), lower = 0, upper = chromosome_length)$val

  E_I <- 10^ibd_log10_pr_cpp(ibd_state_by_v = i$ibd_state_by_v,
                          ibd_by_locus = ibd_state,
                          recombination_rate_by_locus = numeric(),
                          number_of_transmissions = i$number_of_relevant_transmissions,
                          fixed_transmission_masks = i$relevant_masks,
                          pr_v_constant = pr_v_constant,
                          pr_v_biased = pr_v)

  scale <- if(fraction) 1.0 else chromosome_length

  mean <- E_I * scale
  var <- scale^2 * E_I_squared - (E_I*scale)^2
  sd <- sqrt(var)

  list(mean = mean,
       variance = var,
       sd = sd)
}
