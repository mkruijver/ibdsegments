test_that("compare to manual values for ibd", {

  expect_equal(get_ibd_state_2p(c(1,1,1,1), 1, 2), 2)
  expect_equal(get_ibd_state_2p(c(1,1,1,2), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,1,2,1), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,1,2,2), 1, 2), 0)

  expect_equal(get_ibd_state_2p(c(1,2,1,1), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,1,2), 1, 2), 2)
  expect_equal(get_ibd_state_2p(c(1,2,1,3), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,2,1), 1, 2), 2)
  expect_equal(get_ibd_state_2p(c(1,2,2,2), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,2,3), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,3,1), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,3,2), 1, 2), 1)
  expect_equal(get_ibd_state_2p(c(1,2,3,3), 1, 2), 0)

  expect_equal(get_ibd_state_2p(c(1,2,3,4), 1, 2), 0)

})


test_that("compare to manual values for kappa", {

  expect_equal(get_kappa_state(c(1,1,1,1), 1, 2), 0)
  expect_equal(get_kappa_state(c(1,1,1,2), 1, 2), 0)
  expect_equal(get_kappa_state(c(1,1,2,1), 1, 2), 0)
  expect_equal(get_kappa_state(c(1,1,2,2), 1, 2), 0)

  expect_equal(get_kappa_state(c(1,2,1,1), 1, 2), 0)
  expect_equal(get_kappa_state(c(1,2,1,2), 1, 2), 2)
  expect_equal(get_kappa_state(c(1,2,1,3), 1, 2), 1)
  expect_equal(get_kappa_state(c(1,2,2,1), 1, 2), 2)
  expect_equal(get_kappa_state(c(1,2,2,2), 1, 2), 0)
  expect_equal(get_kappa_state(c(1,2,2,3), 1, 2), 1)
  expect_equal(get_kappa_state(c(1,2,3,1), 1, 2), 1)
  expect_equal(get_kappa_state(c(1,2,3,2), 1, 2), 1)
  expect_equal(get_kappa_state(c(1,2,3,3), 1, 2), 0)

  expect_equal(get_kappa_state(c(1,2,3,4), 1, 2), 0)

})
test_that("compare to manual values for multi-person ibd", {

  expect_equal(get_joint_ibd_state(c(1,1,1,1,1,1), 1:3), 2)

  expect_equal(get_joint_ibd_state(c(1,1,1,1,1,2), 1:3), 1)
  expect_equal(get_joint_ibd_state(c(1,1,1,2,1,1), 1:3), 1)
  expect_equal(get_joint_ibd_state(c(1,2,1,1,1,1), 1:3), 1)

  expect_equal(get_joint_ibd_state(c(1,2,1,2,1,2), 1:3), 2)
  expect_equal(get_joint_ibd_state(c(1,2,1,2,1,2), 1:3), 2)
  expect_equal(get_joint_ibd_state(c(1,2,1,2,1,2), 1:3), 2)

  expect_equal(get_joint_ibd_state(c(1,2,1,1,1,2), 1:3), 1)
  expect_equal(get_joint_ibd_state(c(1,2,1,1,1,2), 1:3), 1)
  expect_equal(get_joint_ibd_state(c(1,2,1,1,1,2), 1:3), 1)

  expect_equal(get_joint_ibd_state(c(1,2,1,1,2,3), 1:3), 0)
  expect_equal(get_joint_ibd_state(c(1,2,1,2,1,1), 1:3), 1)
  expect_equal(get_joint_ibd_state(c(1,2,1,2,2,2), 1:3), 1)

  expect_equal(get_joint_ibd_state(c(1,1,2,2,1,1), 1:3), 0)
  expect_equal(get_joint_ibd_state(c(1,1,2,2,2,2), 1:3), 0)
})


test_that("compare to manual values for Jacquard states", {

  # 1: all the same
  expect_equal(get_Jacquard_state(c(1,1,1,1), 1, 2), 1)

  # 2: a=b, c=d
  expect_equal(get_Jacquard_state(c(1,1,2,2), 1, 2), 2)

  # 3: a=b=c or a=b=d
  expect_equal(get_Jacquard_state(c(1,1,1,2), 1, 2), 3)
  expect_equal(get_Jacquard_state(c(1,1,2,1), 1, 2), 3)

  # 4: a=b
  expect_equal(get_Jacquard_state(c(1,1,2,3), 1, 2), 4)

  # 5: a=c=d or b=c=d
  expect_equal(get_Jacquard_state(c(1,2,1,1), 1, 2), 5)
  expect_equal(get_Jacquard_state(c(1,2,2,2), 1, 2), 5)

  # 6: c=d
  expect_equal(get_Jacquard_state(c(1,2,3,3), 1, 2), 6)

  # 7: a=c, b=d or a=d, b=c
  expect_equal(get_Jacquard_state(c(1,2,1,2), 1, 2), 7)
  expect_equal(get_Jacquard_state(c(1,2,2,1), 1, 2), 7)

  # 8: a=c or a=d or b=c or b=d
  expect_equal(get_Jacquard_state(c(1,2,1,3), 1, 2), 8)
  expect_equal(get_Jacquard_state(c(1,2,3,1), 1, 2), 8)
  expect_equal(get_Jacquard_state(c(1,2,2,3), 1, 2), 8)
  expect_equal(get_Jacquard_state(c(1,2,3,2), 1, 2), 8)

  # 9: all different
  expect_equal(get_Jacquard_state(c(1,2,3,4), 1, 2), 9)

})
