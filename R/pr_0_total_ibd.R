#' Probability of no IBD across one or more chromosomes
#'
#' Donnelly (1983) studies the probability that relatives of certain types have
#' no segments of a chromosome in common and provides expressions that can be
#' efficiently computed.
#'
#' @details
#' Types of relationships supported are:
#' \describe{
#' \item{\code{cousin}}{Use `degree = a` and `removal = b` for `a`'th cousins `b` times removed with degree at least one. Default is `degree = 1`.}
#' \item{\code{halfsib}}{Use `removal = 0` (default) for half-siblings and `removal = a` for half-cousins `a` times removed.}
#' \item{\code{grandparent}}{`degree = 1` (default) for grandparents, 2 for great-grandparents and so on.}
#' \item{\code{uncle}}{Use `degree = 0` (default) for uncle and `degree = d` for great^d uncle}
#' }
#'
#' @param relationship_type One of \code{"cousin"}, \code{"halfsib"},
#'                          \code{"grandparent"} or \code{"uncle"}.
#' @param degree See details.
#' @param removal See details.
#' @param removal1 See details.
#' @param removal2 See details.
#' @param chromosome_length Default is 267.77 cM (an estimate of the length of chromosome 1).
#' @param log10 Should the log10 probability be returned? Default is `FALSE`.
#' @return Numeric
#' @references Donnelly K. P. (1983). The probability that related individuals
#'             share some section of genome identical by descent.
#'             Theoretical population biology, 23(1), 34–63.
#'             https://doi.org/10.1016/0040-5809(83)90004-7
#'
#' @seealso [pr_all_genes_survive] for the probability that all autosomal genes
#'          are passed on to the next generation (offspring column in Table 1
#'          of Donnelly (1983))
#' @examples
#' ## Cousin-type: third cousins on a chromosome of length 100 cM
#' degree <- 3
#'
#' p0_3C <- pr_0_total_ibd("cousin", degree = 3, chromosome_length = 100)
#' p0_3C
#'
#' # verify
#' p0_3C_manual <- d_cibd(x = 100, ibd = 0,
#'                        pedigree = pedtools::cousinPed(degree = 3))
#'
#' p0_3C_manual
#' stopifnot(all.equal(p0_3C, p0_3C_manual))
#'
#' ## Half-sib type: half-cousins twice removed
#' p0_H1C_2R <- pr_0_total_ibd("halfsib",
#'                             degree = 1, removal = 2, chromosome_length = 100)
#' p0_H1C_2R
#'
#' p0_H1C_2R_manual <- d_cibd(x = 100, ibd = 0,
#'                            pedigree = pedtools::halfCousinPed(removal = 2))
#'
#' p0_H1C_2R_manual
#'
#' stopifnot(all.equal(p0_H1C_2R, p0_H1C_2R_manual))
#'
#' ## Grandparent-type: great grandparents (degree = 2)
#' p0_GGP <- pr_0_total_ibd("grandparent", degree = 2, chromosome_length = 100)
#' p0_GGP
#'
#' # GGP is a third generation ancestor so n = 3
#' p0_GGP_manual <- d_cibd(x = 100, ibd = 0,
#'                     pedigree = pedtools::linearPed(n = 3),
#'                     ids = c(1, pedtools::leaves(pedtools::linearPed(n = 3))))
#'
#' p0_GGP_manual
#'
#' stopifnot(all.equal(p0_GGP, p0_GGP_manual))
#'
#' ## Uncle-type: degree = 0 for uncle
#' p0_AV <- pr_0_total_ibd("uncle", chromosome_length = 100)
#' p0_AV
#'
#' p0_AV_manual <- d_cibd(x = 100, ibd = 0, pedigree = pedtools::avuncularPed())
#'
#' p0_AV_manual
#'
#' stopifnot(all.equal(p0_AV, p0_AV_manual))
#'
#' ## Reproduce Table 1 of Donnelly (1983)
#' # (historic) chromosome lengths (cM) used in Donnelly (1983)
#' L <- 33 * c(9.12, 8.53, 7.16, 6.59, 6.15, 5.87, 5.31, 4.92, 4.81, 4.71, 4.60,
#'             4.47, 3.56, 3.60, 3.40, 3.20, 3.12, 2.72, 2.48, 2.27, 1.77, 1.64)
#' k <- 4:15
#'
#' tab1 <- data.frame(k=k)
#'
#' tab1$cousin <- pr_0_total_ibd(relationship_type = "cousin",
#'                               degree = rep(2:7, each = 2),
#'                               removal = rep(0:1, times = 6),
#'                               chromosome_length = L)
#'
#' tab1$uncle <- pr_0_total_ibd(relationship_type = "uncle",
#'                              degree = k - 1, chromosome_length = L)
#'
#' # Note the removal on one side only
#' tab1$halfsib <- pr_0_total_ibd(relationship_type = "halfsib",
#'                                removal1 = k - 1,
#'                                removal2 = rep(0, length(k)),
#'                                chromosome_length = L)
#'
#' # Poisson approximation
#' tab1$`exp(-33 * k / 2^k)` <- exp(-33 * k / 2^k)
#'
#' # Note that k corresponds to great^k grandparent,
#' # i.e. a (k+2)'th generation ancestor
#' # (not great^(k-1) and (k+1)'th generation ancestor as printed)
#'
#' tab1$grandparent <- pr_0_total_ibd(relationship_type = "grandparent",
#'                                    degree = k, chromosome_length = L)
#'
#' tab1
#' @export
pr_0_total_ibd <- function(relationship_type = c("cousin", "grandparent",
                                                 "halfsib", "uncle"),
                      degree,
                      removal,
                      removal1,
                      removal2,
                      chromosome_length,
                      log10 = FALSE)
{

  if (relationship_type == "grandparent"){
    if (missing(degree)) degree <- 1L
    .validate_integer(degree, "degree")

    if (!missing(removal) || !missing(removal1) || !missing(removal2)){
      stop("Only degree argument (not removals) is used for grandparent type")
    }
    if (any(degree < 1)) stop("degree should be at least one")

    # vectorise on degree parameter
    if (length(degree) > 1){
      return(sapply(degree, function(degree)
        pr_0_total_ibd(relationship_type = relationship_type, degree = degree,
                       chromosome_length = chromosome_length, log10 = log10)))
    }

    ctmc <- ctmc_donnelly_grandparent(d = degree)
  }
  else if (relationship_type == "cousin"){
    if (missing(degree)) degree <- 1L
    .validate_integer(degree, "degree")
    if (missing(removal)) removal <- 0L

    .validate_integer(removal, "removal")

    if (any(degree < 1)) stop("degree should be at least one")
    if (any(removal < 0)) stop("removal should not be negative")

    if (!missing(removal1) || !missing(removal2)){
      stop("Only degree and removal arguments (not removal1) is used for cousin type")
    }

    # vectorise over degree and removal
    if (length(degree)>1){
      if (length(degree) != length(removal)){
        stop("degree and removal need to have the same length")
      }

      return(sapply(seq_along(degree), function(i){
        pr_0_total_ibd(relationship_type = relationship_type, degree = degree[i],
                       removal = removal[i],
                       chromosome_length = chromosome_length, log10 = log10)
      }))
    }

    ctmc <- ctmc_donnelly_cousin(degree = degree, removal = removal)
  }
  else if (relationship_type == "halfsib"){
    if (!missing(removal) && (!missing(removal1) || !missing(removal2))){
      stop("Either provide removal or removal1 but not both")
    }
    if (missing(removal)) removal <- 0L

    if (missing(removal1)) removal1 <- removal
    if (missing(removal2)) removal2 <- removal

    # vectorise over removal1 and removal2
    if (length(removal1)>1){
      if (length(removal1) != length(removal2)){
        stop("removal1 and removal2 need to have the same length")
      }

      return(sapply(seq_along(removal1), function(i){
        pr_0_total_ibd(relationship_type = relationship_type,
                       removal1 = removal1[i],
                       removal2 = removal2[i],
                       chromosome_length = chromosome_length, log10 = log10)
      }))
    }

    if (removal1 < 0) stop("removal1 should not be negative")
    if (removal2 < 0) stop("removal2 should not be negative")

    ctmc <- ctmc_donnelly_halfsib(removal1, removal2)
  }
  else if (relationship_type == "uncle"){
    if (missing(degree)) degree <- 0L
    .validate_integer(degree, "degree")

    if (!missing(removal) || !missing(removal1) || !missing(removal2)){
      stop("Only degree argument (not removals) is used for uncle type")
    }

    # vectorise on degree parameter
    if (length(degree) > 1){
      return(sapply(degree, function(degree)
        pr_0_total_ibd(relationship_type = relationship_type, degree = degree,
                       chromosome_length = chromosome_length, log10 = log10)))
    }

    if (degree < 0) stop("degree should not be negative")

    ctmc <- ctmc_donnelly_uncle(degree)
  }
  else{
    stop("Relationship ", relationship_type, " is not supported")
  }

  log10p <- sum(sapply(chromosome_length, function(l){
    log10(sum(expm::expm.Higham08(l/100 * ctmc$Q_dishonest) %*% ctmc$pi0_dishonest))
  }))

  if (log10) log10p else 10 ^ log10p
}

ctmc_donnelly_grandparent <- function(d){
  # d = 1 means grandparent (2nd generation ancestor)
  # In Table 1, k = n, (not k=n+1 as printed!)
  # so k=4 means n=4 aka the 6th gen. ancestor (n+2)'th
  # so d = k

  Q_diag <- rep(-d, d+1)
  Q_super_diag <- 1:d
  Q_sub_diag <- d:1

  Q <- matrix(0, d + 1, d + 1)
  diag(Q) <- Q_diag

  Q[row(Q) == col(Q) - 1] <- Q_super_diag
  Q[row(Q) == col(Q) + 1] <- Q_sub_diag

  pi0 <- choose(d, 0:d) / (2^d)

  # first orbit is form hitting set
  Q_dishonest = Q[-1, -1, drop = FALSE]
  pi0_dishonest <- pi0[-1, drop = FALSE]

  return(list(Q = Q, pi0 = pi0,
              Q_dishonest = Q_dishonest,
              pi0_dishonest = pi0_dishonest))
}

ctmc_donnelly_cousin <- function(degree, removal){

  d <- 2 * degree + removal + 4

  M <- matrix(c(0,2,2,0,2,0,0,
                2,0,0,2,0,1,0,
                2,0,0,2,0,1,0,
                0,2,2,0,0,0,2,
                2,0,0,0,0,2,0,
                0,2,2,0,4,0,4,
                0,0,0,2,0,2,0), nrow = 7, ncol = 7, byrow = TRUE)

  Q0 <- M - d * diag(7)

  number_of_blocks <- d - 6 + 1

  Q <- ctmc_donnelly_block_Q(M = M, d = d, number_of_blocks = number_of_blocks,
                             first_sub = -6, first_super = 1)

  # u is of length 7
  u <- c(rep(1,5), 2, 1) / 8
  # initial vector has length d - 6 + 1 = 5

  u_rep <- rep(u, d - 6 + 1)
  binom_rep <- rep(choose(d-6, 0:(d-6)), each = length(u))

  pi0 <- binom_rep * u_rep / (2^(d-6))

  Q_dishonest <- Q[-(1:2),-(1:2)]
  pi0_dishonest <- pi0[-(1:2)]

  list(Q = Q,
       pi0 = pi0,
       Q_dishonest = Q_dishonest,
       pi0_dishonest = pi0_dishonest)
}

ctmc_donnelly_halfsib <- function(removal1, removal2){

  d <- removal1 + removal2 + 2

  M <- matrix(c(0, 2,
                2, 0),
              nrow = 2, ncol = 2, byrow = TRUE)

  number_of_blocks <- d - 1
  Q <- ctmc_donnelly_block_Q(M = M, d = d, number_of_blocks = number_of_blocks,
                             first_sub = -2, first_super = 1)

  # u is of length 2
  u <- 1/4 * c(2, 2)

  u_rep <- rep(u, d - 2 + 1)
  binom_rep <- rep(choose(d - 2, 0 : (d-2)), each = length(u))

  pi0 <- binom_rep * u_rep / (2^(d-2))

  # first orbit is hitting set
  Q_dishonest <- Q[-1,-1, drop = FALSE]
  pi0_dishonest <- pi0[-1, drop = FALSE]

  list(Q = Q,
       pi0 = pi0,
       Q_dishonest = Q_dishonest,
       pi0_dishonest = pi0_dishonest)
}

ctmc_donnelly_uncle <- function(degree = 0){

  d <- 5 + degree

  M <- matrix(c(1, 2, 0, 2,
                2, 0, 2, 1,
                0, 2, 1, 2,
                2, 1, 2, 0),
              nrow = 4, ncol = 4, byrow = TRUE)

  number_of_blocks <- d - 5 + 1
  Q <- ctmc_donnelly_block_Q(M = M, d = d, number_of_blocks = number_of_blocks,
                             first_sub = -5, first_super = 1)

  # u is of length 4
  u <- 1/4 * c(1, 1, 1, 1)

  u_rep <- rep(u, d - 5 + 1)
  binom_rep <- rep(choose(d - 5, 0 : (d-5)), each = length(u))

  pi0 <- binom_rep * u_rep / (2^(d-5))

  # first orbit is hitting set
  Q_dishonest <- Q[-(1:2),-(1:2), drop = FALSE]
  pi0_dishonest <- pi0[-(1:2), drop = FALSE]

  list(Q = Q,
       pi0 = pi0,
       Q_dishonest = Q_dishonest,
       pi0_dishonest = pi0_dishonest)
}

ctmc_donnelly_block_Q <- function(M, d, number_of_blocks,
                                  first_sub = -2, first_super = 1){
  # first_super  = -1 means that first block above diagonal is equals 1 * I
  # first_sub = -2 means that the first sub diagonal block is equals (d - 2 ) * I
  # e.g. block 2 (below block 1) is (d - 2) I

  m <- nrow(M)
  I <- diag(m)

  Q <- matrix(data = 0,
              nrow = m * number_of_blocks,
              ncol = m * number_of_blocks)

  Q0 <- M - d * I

  # block 1
  Q[1:m, 1:m] <- Q0

  if (number_of_blocks >= 2){
    for (i in 2:(number_of_blocks)){
      idx <- (1 + (i-1) * m) : (i*m)

      Q[idx, idx] <- Q0
      Q[idx - m, idx] <- (first_super + (i - 2)) * I

      Q[idx, idx - m] <- (d + first_sub - (i - 2)) * I
    }
  }

  Q
}

half_sib_type_ped <- function(removal1 = 0, removal2 = 0){
  # half sib type pedigree following Donnelly (1983)
  # in pedtools this would be called a half cousin pedigree
  # pedtools::halfCousinPed does not support removal1 != removal2
  x = pedtools::linearPed(1 + removal1, sex = 1)
  y = pedtools::linearPed(1 + removal2)
  pedtools::mergePed(x, y, by = 1, relabel = TRUE)
}
