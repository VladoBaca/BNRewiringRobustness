context("Utility functions used around the whole package")

test_that("all_binary_vectors returns correct vectors", {
  expect_equal(all_binary_vectors(0), list(numeric(0)))
  expect_equal(all_binary_vectors(1), list(c(0), c(1)))
  expect_equivalent(all_binary_vectors(2), list(c(0,0), c(0,1), c(1,0), c(1,1)))
})

test_that("binary_vector_to_number and binary_vector_to_one_based_index return correct numbers", {
  expect_equal(binary_vector_to_number(c(0,1,1,0,1)), 13)
  expect_equal(binary_vector_to_one_based_index(c(0,0,1,0,1)), 6)
})

test_that("pre_to_post_vector return correct vector", {
  expect_equal(pre_to_post_vector(all_binary_vectors(3)), c(0,4,2,6,1,5,3,7))
})
