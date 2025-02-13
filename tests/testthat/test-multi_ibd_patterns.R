

test_that("verify against ribd", {
  test_peds <- list(pedtools::nuclearPed(nch = 2),
                      pedtools::avuncularPed(),
                      pedtools::cousinPed(degree = 1),
                      pedtools::cousinPed(degree = 2),
                      pedtools::cousinPed(degree = 5))

  for (ped in test_peds){
    expected <- ribd::multiPersonIBD(ped,
                                     ids = pedtools::leaves(ped))
    observed <- multi_ibd_patterns(ped)

    expect_equal(observed, expected, ignore_attr = TRUE)
  }
})
