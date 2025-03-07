test_that("total_ibd_dist_moments agrees with total_ibd_dist integrated", {


  L1 <- 123
  L2 <- 231

  pedigrees <- list(pedtools::cousinPed(degree = 3),
                    pedtools::fullSibMating(1))

  for (pedigree in pedigrees){

    # one chromosome
    dist1 <- total_ibd_dist(pedigree = pedigree, fraction = FALSE,
                            ibd_state = 1, chromosome_length = L1)
    moments1 <- total_ibd_dist_moments(pedigree, chromosome_length = L1)

    mean1 <- E(dist1)
    var1 <- var(dist1)
    sd1 <- sd(dist1)

    expect_equal(moments1$mean, mean1)
    expect_equal(moments1$var, var1)
    expect_equal(moments1$sd, sd1)

    # second chromosome
    dist2 <- total_ibd_dist(pedigree = pedigree, fraction = FALSE,
                            ibd_state = 1, chromosome_length = L2)

    moments2 <- total_ibd_dist_moments(pedigree, chromosome_length = L2)

    mean2 <- E(dist2)
    var2 <- var(dist2)
    sd2 <- sd(dist2)

    expect_equal(moments2$mean, mean2)
    expect_equal(moments2$var, var2)
    expect_equal(moments2$sd, sd2)

    # convolution
    moments12 <- total_ibd_dist_moments(pedigree, chromosome_length = c(L1, L2))

    expect_equal(moments12$mean, mean1 + mean2)
    expect_equal(moments12$variance, var1 + var2)

    # fraction
    moments_fraction12 <- total_ibd_dist_moments(pedigree,
                                                 chromosome_length = c(L1, L2),
                                                 fraction = TRUE)

    mean_fraction12 <- (mean1 + mean2) / (L1 + L2)
    var_fraction12 <- (var1 + var2) / (L1 + L2)^2

    expect_equal(moments_fraction12$mean, mean_fraction12)
    expect_equal(moments_fraction12$var, var_fraction12)
    expect_equal(moments_fraction12$sd, sqrt(var_fraction12))
  }

})
