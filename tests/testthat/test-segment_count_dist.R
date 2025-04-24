
test_that("verify segment_count_dist pr(0) against pr_0_total_ibd", {

  # cousin pedigrees
  for (degree in 1:4){
    expected <- pr_0_total_ibd(relationship_type = "cousin", degree = degree,
                   chromosome_length = 123)

    dist <- segment_count_dist(pedigree = pedtools::cousinPed(degree = degree),
                       chromosome_length = 123)

    observed <- unname(dist$px[match(0, dist$x)])


    expect_equal(observed, expected)
  }

})

test_that("verify segment_count_dist moments against simulations", {

  # cousin pedigrees
  for (degree in 1:4){

    pedigree <- pedtools::cousinPed(degree = degree)

    r <- r_cibd(pedigree = pedigree, chromosome_length = 123, n = 5e3)

    sd(r$stats$segment_count)

    dist <- segment_count_dist(pedigree = pedigree,
                               chromosome_length = 123)

    E_dist <- E(dist)
    expect_equal(mean(r$stats$segment_count), E_dist, tolerance = .25)

    sd_dist <- sd(dist)
    expect_equal(sd(r$stats$segment_count), sd_dist, tolerance = .25)
  }

})
