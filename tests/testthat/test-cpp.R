test_that("C++ unit tests pass", {
  result <- skindiff:::.cpp_run_tests()
  expect_true(isTRUE(as.logical(result)))
})
