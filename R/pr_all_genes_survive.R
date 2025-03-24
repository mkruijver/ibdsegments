#' Probability of passing down all autosomal genes to offspring
#'
#' Donnelly (1983) presents a way to efficiently calculate the probability of
#' passing down all genes to the next generation. This includes both the
#' paternally and maternally inherited genes, so at least two offspring are
#' needed for this probability to be non-zero.
#'
#' @param number_of_children How many offspring? Result is vectorised over this
#'                           parameter.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of
#'                          chromosome 1). `length(chromosome_length)>1`, the
#'                          probability of passing down all chromosomes will be
#'                          computed.
#' @param log10 Should the log10 probability be returned? Default is `FALSE`.
#' @return Numeric
#' @references Donnelly K. P. (1983). The probability that related individuals
#'             share some section of genome identical by descent.
#'             Theoretical population biology, 23(1), 34–63.
#'             https://doi.org/10.1016/0040-5809(83)90004-7
#' @examples
#' ## Example from Donnelly (1983)
#' # (historic) chromosome lengths (cM) used in Donnelly (1983)
#' L <- 33 * c(9.12, 8.53, 7.16, 6.59, 6.15, 5.87, 5.31, 4.92, 4.81, 4.71, 4.60,
#'             4.47, 3.56, 3.60, 3.40, 3.20, 3.12, 2.72, 2.48, 2.27, 1.77, 1.64)
#'
#' # Reproduce the "Offspring" column in Table 1
#' number_of_children <- 1:16
#' pr_survival <- pr_all_genes_survive(number_of_children, chromosome_length = L)
#'
#' # Plot the results (compare to Fig. 10)
#' plot(number_of_children, pr_survival)
#'
#' # Add the Poisson approximation
#' k <- 2:17 - 1
#' lines(k+1, exp(-k*33 / (2^k)), lty=2)
#'
#' ## Example using general purpose methods
#' # The probability of passing down all genes to k offspring is equal to the
#' # probability that the joint IBD state of k half siblings is 0 everywhere -
#' # there is no point where they all inherited the same DNA from the common
#' # parent.
#'
#' # Define a function to create a pedigree of half siblings
#' ped_halfsibs <- function(number_of_half_sibs){
#'   ped <- pedtools::nuclearPed(nch = 1)
#'   for (k in 2:number_of_half_sibs) {
#'     ped <- pedtools::addChild(ped, parents = c(1), verbose = FALSE)
#'   }
#'   ped
#' }
#'
#' # Compute the probability that a chromosome of length 100 cm survives
#' # if the next generation consists of 5 children
#' p1 <- pr_all_genes_survive(number_of_children = 5L, chromosome_length = 100)
#' p2 <- d_cibd(x = 100, ibd = 0, pedigree = ped_halfsibs(5))
#'
#' p1
#' p2
#' stopifnot(all.equal(p1, p2))
#'
#' # If you have five children, which fraction of chromosome 1
#' # is passed down to the next generation? Assume a length of 268 cM
#'
#' # Pr. of passing down the complete grand-maternal and grand-paternal parts
#' # is 42%
#' pr_all_genes_survive(number_of_children = 5L, chromosome_length = 268)
#'
#' # The distribution of the fraction that is passed down has a point
#' # mass of 0.42 at 100% and has a continuous density with weight 0.58
#' dist <- total_ibd_dist(ped_halfsibs(5), ibd_state = 0,
#' chromosome_length = 267, fraction = TRUE)
#' dist
#'
#' plot(dist)
#' @export
pr_all_genes_survive <- function(number_of_children, chromosome_length = 267.66,
                                 log10 = FALSE){

  .validate_integer(number_of_children, "number_of_children")
  .validate_chromosome_length(chromosome_length)
  .validate_logical(log10, "log10")

  if (length(number_of_children) > 1){
    return(sapply(number_of_children, function(n)
      pr_all_genes_survive(n, chromosome_length = chromosome_length,
                           log10 = log10)))
  }

  if (number_of_children < 2){
    return(if(log10) log10(0.) else 0.)
  }

  ctmc <- ctmc_offspring(number_of_children)

  log10p <- sum(sapply(chromosome_length, function(l){
    log10(sum(expm::expm.Higham08(l/100 * ctmc$Q_dishonest) %*% ctmc$pi0_dishonest))
  }))

  if (log10) log10p else 10 ^ log10p
}

ctmc_offspring <- function(d){

  # Donnelly (1983) actually uses only d/2+1 orbits if d is even
  # and (d+1)/2 orbits if d is odd

  # instead we simply use a state space that
  # represents the number of 1's like the grandparent case
  # then the hitting set is 0...0 and 1...1

  number_of_orbits <- d + 1

  Q_diag <- rep(-d, number_of_orbits)
  Q_super_diag <- 1:d
  Q_sub_diag <- d:1

  Q <- matrix(0, number_of_orbits, number_of_orbits)
  diag(Q) <- Q_diag

  Q[row(Q) == col(Q) - 1] <- Q_super_diag
  Q[row(Q) == col(Q) + 1] <- Q_sub_diag

  pi0 <- choose(d, 0:d) / (2^d)

  Q_dishonest = Q[-c(1,d+1), -c(1,d+1), drop = FALSE]
  pi0_dishonest <- pi0[-c(1,d+1), drop = FALSE]

  return(list(Q = Q, pi0 = pi0,
              Q_dishonest = Q_dishonest,
              pi0_dishonest = pi0_dishonest))
}
