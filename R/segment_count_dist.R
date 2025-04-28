#' Compute probability distribution of IBD segment count
#'
#' The `segment_count_dist` function computes the probability distribution
#' of the number of IBD segments
#' (i.e., the count of IBD intervals) between pairs of individuals
#' in a pedigree, for a given IBD state and chromosome length(s).
#'
#' This function is analogous to [total_ibd_dist()] but focuses on
#' the number of IBD segments rather than the total IBD length.
#'
#' @param pedigree Pedigree in [`pedtools::ped`] form.
#' @param ids Ids for which IBD is observed. Default is `pedtools::leaves(pedigree)`.
#' @param states One of `"ibd"` (default), `"kappa"`, `"identity"` or `"detailed"`.
#' @param ibd_state Default is 1.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#' @param convolve Should the distribution of the sum (across chromosomes) be obtained?
#' @return object of class `segment_count_dist`
#'
#' @seealso [total_ibd_dist()] for the distribution of IBD *length*.
#'
#' @examples
#' ped_hs <- pedtools::halfSibPed()
#'
#' dist <- segment_count_dist(ped_hs, chromosome_length = 100)
#' dist
#' plot(dist)
#' E(dist)
#' sd(dist)
#'
#' r <- r_cibd(n = 1e4, pedigree = ped_hs, chromosome_length = 100)
#' mean(r$stats$segment_count)
#' sd(r$stats$segment_count)
#'
#' @importFrom stats qpois dpois
#' @export
segment_count_dist <- function(pedigree,
                           ids = pedtools::leaves(pedigree),
                           states = "ibd",
                           ibd_state = 1L,
                           chromosome_length = 267.77,
                           convolve = TRUE){

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

  # compute the probability distribution of k = 0, 1, 2, ... n_max+1 IBD segments
  # if there are n = 0, 1, 2, ..., n_max recombinations
  V <- pr_number_of_segments_by_n(ibd_state = ibd_state,
                                            ibd_state_by_v = i$ibd_state_by_v, n_max = joint_n_max,
                                            number_of_transmissions = i$number_of_relevant_transmissions,
                                            masks = i$relevant_masks)
  # colnames(V) <- paste0("n=", 0:(ncol(V) - 1))
  # rownames(V) <- paste0("k=", 0:(nrow(V) - 1))

  chromosome_dists <- lapply(seq_along(chromosome_length),
                             function(i_chromosome){
    lambda <- 0.01 * chromosome_length[i_chromosome] *
      i$number_of_relevant_transmissions

    n_max <- stats::qpois(1.0 - 1e-16, lambda = lambda)
    n <- 0:n_max

    pr_n <- stats::dpois(x = n, lambda = lambda)
    pr_n_padded <- c(pr_n, rep(0, ncol(V) - length(pr_n)))

    x <- 0:(nrow(V) - 1)
    px <- colSums(t(V) * pr_n_padded)

    idx_pos <- px > 0

    dist_chromosome <- list(x = x[idx_pos],
                            px = px[idx_pos])

    class(dist_chromosome) <- "segment_count_dist"

    dist_chromosome
  })

  if (length(chromosome_dists) > 1 && convolve){
    return(convolve_segment_count_dists(chromosome_dists))
  }

  if (length(chromosome_length) == 1){
    chromosome_dists <- chromosome_dists[[1]]
  }

  chromosome_dists
}

#' @method E segment_count_dist
#' @export
E.segment_count_dist <- function(x, m = 1, ...){
  sum(x$x^m * x$px)
}

#' @method var segment_count_dist
#' @export
var.segment_count_dist <- function(x, ...){
  E.segment_count_dist(x, m = 2) - E.segment_count_dist(x)^2
}

#' @method sd segment_count_dist
#' @export
sd.segment_count_dist <- function(x, ...){
  sqrt(var.segment_count_dist(x))
}

#' @method d segment_count_dist
#' @export
d.segment_count_dist <- function(x, ...){

  dist <- x

  function(x){

    idx <- match(x, dist$x)

    pr <- dist$px[idx]
    pr[is.na(pr)] <- 0.

    pr
  }
}


#' @importFrom graphics plot points
#' @method plot segment_count_dist
#' @export
plot.segment_count_dist <- function(x, ...){

  graphics::plot(x$x, x$px, type="h",
       xlab = "Segment count",
       ylab = "Probability", ...)
  graphics::points(x$x, x$px)
}
