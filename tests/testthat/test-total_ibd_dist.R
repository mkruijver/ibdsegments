
test_that("total_ibd_dist point masses agree with d_cibd", {


  L <- 123

  pedigrees <- list(pedtools::cousinPed(degree = 3),
                    pedtools::fullSibMating(1))

  for (pedigree in pedigrees){
    d0 <- total_ibd_dist(pedigree = pedigree, fraction = FALSE,
                        ibd_state = 0, chromosome_length = L)
    p <- d0$point_mass$px[2]
    p_cibd <-  d_cibd(x = L, ibd = 0, pedigree = pedigree)

    expect_equal(p, p_cibd)

    d1 <- total_ibd_dist(pedigree = pedigree, fraction = FALSE,
                         ibd_state = 1, chromosome_length = L)

    p1 <- d1$point_mass$px[2]
    p1_cibd <- d_cibd(x = L, ibd = 1, pedigree = pedigree)

    expect_equal(p1, p1_cibd)
  }
})

test_that("total_ibd_dist agrees with Ball and Ivanov", {

  # Reproduce part of Table 3 in Ball and Ivanov (2005)
  # doi: 10.1016/j.mbs.2005.04.005

  expected <- matrix(c(
      0.0677, 0.0682, 0.0732, 0.0968, 0.1298, 0.1666, 0.2070,
      0.2760, 0.2773, 0.2897, 0.3442, 0.4106, 0.4744, 0.5350,
      0.5095, 0.5110, 0.5250, 0.5841, 0.6497, 0.7070, 0.7567,
      0.6920, 0.6932, 0.7052, 0.7538, 0.8045, 0.8456, 0.8789),
    byrow = TRUE, nrow = 4)

  m_rows <- 2:5
  s_cols <- c(0.01, 0.1, 1, 5, 10, 15, 20)

  # make half sib pedigrees
  peds_hs <- list()
  peds_hs[[2]] <- pedtools::halfSibPed()
  for(m in m_rows[-1]){
    peds_hs[[m]] <- pedtools::addChildren(peds_hs[[m - 1]],
                                          father = 2, verbose = FALSE)
  }

  # fill matrix
  observed <- matrix(nrow = nrow(expected), ncol = ncol(expected))
  for (i_row in seq_len(nrow(expected))){

    m <- m_rows[i_row]

    dist <- total_ibd_dist(pedigree = peds_hs[[m]], chromosome_length = 100)
    f <- d(dist)
    w <- dist$weight_continuous
    p <- function(x) dist$point_mass$px[1] +
                     w*integrate(f = f, lower = 0, upper = x)$val

    for (i_col in seq_len(ncol(expected))){
      s <- s_cols[i_col]

      observed[i_row, i_col] <- p(s)
    }
  }

  testthat::expect_equal(observed, expected, tolerance = 1e-4)
})
