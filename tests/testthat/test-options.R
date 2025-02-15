
test_that("ignore_irrelevant_transmissions option set and get", {

  original_value <- get_option_ignore_irrelevant_transmissions()

  # set to TRUE
  set_option_ignore_irrelevant_transmissions(TRUE)
  expect_equal(get_option_ignore_irrelevant_transmissions(), TRUE)

  # set to FALSE
  set_option_ignore_irrelevant_transmissions(FALSE)
  expect_equal(get_option_ignore_irrelevant_transmissions(), FALSE)

  # restore
  set_option_ignore_irrelevant_transmissions(original_value)
  expect_equal(get_option_ignore_irrelevant_transmissions(), original_value)
})

