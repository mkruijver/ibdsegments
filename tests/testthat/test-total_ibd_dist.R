
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
