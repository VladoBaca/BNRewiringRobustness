context("Functionality for work with PBN instances")

test_that("extract_instance correctly extracts a BN instance", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))

  bn <- extract_instance(pbn, 1)

  expect_equivalent(bn$genes, c("G1", "G2", "G3"))

  expect_equivalent(bn$interactions$G1$input, c(3))
  expect_equivalent(bn$interactions$G1$func, c(0,1))

  expect_equivalent(bn$interactions$G2$input, c(1,2))
  expect_equivalent(bn$interactions$G2$func, c(0,0,1,0))

  expect_equivalent(bn$interactions$G3$input, c(1))
  expect_equivalent(bn$interactions$G3$func, c(0,1))

  expect_equivalent(bn$fixed, c(-1, -1, -1))
  expect_equivalent(class(bn), "BooleanNetwork")
})

test_that("build_transition_table correctly builds the transition table", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  bn <- extract_instance(pbn, 2)

  actual_table <- build_transition_table(bn, c(pbn$genes))

  expected_table <- matrix(c(
    0,1,0,
    1,1,0,
    0,0,0,
    1,0,0,
    0,1,1,
    1,1,1,
    0,1,1,
    1,1,1
  ), nrow = 8, ncol = 3, byrow = TRUE)

  expect_equivalent(actual_table, expected_table)
})
