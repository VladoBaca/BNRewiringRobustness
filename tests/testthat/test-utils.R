context("Utility functions used around the whole package")

test_that("all_binary_vectors returns correct vectors", {
  expect_equal(all_binary_vectors(0), list(numeric(0)))
  expect_equal(all_binary_vectors(1), list(c(0), c(1)))
  expect_equivalent(all_binary_vectors(2), list(c(0,0), c(0,1), c(1,0), c(1,1)))
})
