coefficient_range <- list(kappa = 0:2,
                          identity = 1:9,
                          detailed = 1:15)

coefficient_short_label <- list(kappa = "k",
                                identity = "D",
                                detailed = "d")

compute_single_locus_relatedness_coefficients <- function(pedigree,
                                                          ids = pedtools::leaves(pedigree),
                                                          coefficient){
  values <- coefficient_range[[coefficient]]

  setNames(sapply(values, function(x){
    d_ibd(ibd = x, pedigree = pedigree, ids = ids,
           states = coefficient)
  }),
  nm = paste0(coefficient_short_label[[coefficient]], values))
}

compute_two_locus_relatedness_coefficients <- function(pedigree, recombination_rate,
                                                       ids = pedtools::leaves(pedigree),
                                                       coefficient){


  dim_names <- paste0(coefficient_short_label[[coefficient]],
                      coefficient_range[[coefficient]])
  values <- coefficient_range[[coefficient]]

  probabilities <- matrix(data = NA, nrow = length(values),
                          ncol = length(values),
                          dimnames = list(dim_names, dim_names))

  for (i1 in seq_along(values)){
    for (i2 in seq_along(values)){
      probabilities[i1, i2] <-  d_ibd(c(values[i1], values[i2]),
                                       recombination_rate_by_locus = recombination_rate,
                                       pedigree = pedigree, ids = ids, states = coefficient)
    }
  }

  probabilities
}

common_peds <- list(pedtools::nuclearPed(nch = 2),
                    pedtools::avuncularPed(),
                    pedtools::cousinPed(degree = 1),
                    pedtools::cousinPed(degree = 2),
                    pedtools::cousinPed(degree = 5))

test_that("verify single locus kappa's against ribd", {

  for (founder_inbreeding in c(FALSE, TRUE)){
    for (ped in common_peds){

      if (founder_inbreeding){
        founders <- pedtools::founders(ped)
        ped <- pedtools::setFounderInbreeding(ped, founders,
                                              value = runif(n = length(founders)))
      }

      expected <- ribd::identityCoefs(ped, ids = pedtools::leaves(ped))[9:7]
      observed <- compute_single_locus_relatedness_coefficients(
        pedigree = ped, coefficient = "kappa")

      expect_equal(unname(observed), expected)
    }
  }
})


test_that("verify single locus Jacquard coefficients against ribd", {

  for (founder_inbreeding in c(FALSE, TRUE)){
    for (ped in common_peds){
      if (founder_inbreeding){
        founders <- pedtools::founders(ped)
        ped <- pedtools::setFounderInbreeding(ped, founders,
                                              value = runif(n = length(founders)))
      }

      expected <- ribd::identityCoefs(ped, ids = pedtools::leaves(ped))[1:9]
      observed <- compute_single_locus_relatedness_coefficients(
        pedigree = ped, coefficient = "identity")

      expect_equal(observed, expected, ignore_attr = TRUE)
    }
  }
})


test_that("verify single locus detailed Jacquard coefficients against ribd", {

  for (founder_inbreeding in c(FALSE, TRUE)){
    for (ped in common_peds){
      if (founder_inbreeding){
        founders <- pedtools::founders(ped)
        ped <- pedtools::setFounderInbreeding(ped, founders,
                                              value = runif(n = length(founders)))
      }

      expected <- ribd::identityCoefs(ped, ids = pedtools::leaves(ped), detailed = TRUE)
      observed <- compute_single_locus_relatedness_coefficients(
        pedigree = ped, coefficient = "detailed")

      expect_equal(observed, expected, ignore_attr = TRUE)
    }
  }
})

test_that("verify two locus kappa's against ribd::twoLocusIBD", {

  for (ped in common_peds){
    ids <- pedtools::leaves(ped)

    expected <- ribd::twoLocusIBD(ped, ids = ids, rho = 0.1)
    observed <- compute_two_locus_relatedness_coefficients(ped, recombination_rate = 0.1,
                                                           ids = ids, coefficient = "kappa")

    expect_equal(observed, expected, ignore_attr = TRUE)
  }

})

test_that("verify two locus Jacquard against ribd::twoLocusIdentity", {
  for (ped in common_peds){
    ids <- pedtools::leaves(ped)

    expected <- ribd::twoLocusIdentity(ped, ids = pedtools::leaves(ped),
                                       rho = 0.01)
    observed <- compute_two_locus_relatedness_coefficients(ped,
                                                           recombination_rate = 0.01, ids = ids,
                                                           coefficient = "identity")

    expect_equal(observed, expected)
  }
})

## this is not (yet?) implemented in ribd
# test_that("verify two locus detailed Jacquard against ribd::twoLocusIdentity", {
#   for (ped in common_peds){
#     ids <- pedtools::leaves(ped)
#
#     expected <- ribd::twoLocusIdentity(ped, ids = pedtools::leaves(ped),
#                                        rho = 0.01, detailed = TRUE)
#     observed <- compute_two_locus_relatedness_coefficients(ped,
#                     recombination_rate = 0.01, ids = ids,
#                     states = "detailed")
#
#     expect_equal(observed, expected)
#   }
# })

condensed_by_detailed_states <- list(`1` = 1, `2` = 6, `3` = 2:3,
                                     `4` = 7, `5` = 4:5, `6` = 8, `7` = c(9, 12),
                                     `8` = c(10, 11, 13, 14), `9` = 15)

map_matrix <- t(sapply(condensed_by_detailed_states, function(detailed_states){
  m <- numeric(15)

  m[detailed_states] <- 1
  m
}))

two_locus_detailed_to_condensed <- function(d){
  if ((nrow(d) != 15) | (ncol(d) != 15) | (!is.matrix(d))){
    stop("d needs to be a 15 x 15 matrix")
  }

  map_matrix %*% d %*% t(map_matrix)
}

test_that("verify two locus detailed Jacquard agrees with one locus and collapsed", {

  for (ped in c(common_peds, list(pedtools::fullSibMating(n = 1)))){

    d1 <- compute_single_locus_relatedness_coefficients(pedigree = ped,
                                                        coefficient = "detailed")
    d2 <- compute_two_locus_relatedness_coefficients(pedigree = ped,
                                                     recombination_rate = 0.123, coefficient = "detailed")

    d2_collapsed <- two_locus_detailed_to_condensed(d2)
    d2_condensed <- ribd::twoLocusIdentity(ped, ids = pedtools::leaves(ped),
                                           rho = 0.123)

    expect_equal(d2_collapsed, d2_condensed, ignore_attr = TRUE)
  }

})
