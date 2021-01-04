context("PBN functionality of the package")

test_that("create_pbn_from_rg correctly builds PBN", {
  rg <- load_RG(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  pbn <- create_pbn_from_rg(rg)

  expect_equivalent(pbn$genes, c("G1", "G2", "G3"))
  expect_equivalent(pbn$gene_regulators_indexes, list(c(3), c(1,2), c(1)))
  expect_equivalent(pbn$gene_function_vectors, list(list(c(0, 1)), list(c(0,0,1,0), c(1,0,1,1)), list(c(0, 1))))
  expect_equivalent(pbn$function_index_combinations, list(c(1,1,1), c(1,2,1)))
})

test_that("create_pbn_from_rg correctly validates multiple edges", {
  rg <- load_RG(system.file("examples", "rgs", "test_incorrect_multiple.rg", package = "BNRewiringRobustness"))

  expect_error({create_pbn_from_rg(rg)}, regexp = "Not all the edges are unique in the regulatory graph!")
})

test_that("get_parametrisation_by_index_combination_vector fetches correct parametrisation", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  parametrisation_and <- get_parametrisation_by_index_combination_vector(pbn$gene_function_vectors, c(1,1,1))
  parametrisation_or <- get_parametrisation_by_index_combination_vector(pbn$gene_function_vectors, c(1,2,1))

  expect_equivalent(parametrisation_and, list(c(0, 1), c(0,0,1,0), c(0, 1)))
  expect_equivalent(parametrisation_or, list(c(0, 1), c(1,0,1,1), c(0, 1)))

})

test_that("get_parametrisation_by_index fetches correct parametrisation", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  parametrisation_and <- get_parametrisation_by_index(pbn, 1)
  parametrisation_or <- get_parametrisation_by_index(pbn, 2)

  expect_equivalent(parametrisation_and, list(c(0, 1), c(0,0,1,0), c(0, 1)))
  expect_equivalent(parametrisation_or, list(c(0, 1), c(1,0,1,1), c(0, 1)))
})


test_that("generate_BN correctly generates BN", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  bn <- generate_BN(pbn$genes, pbn$gene_regulators_indexes, get_parametrisation_by_index(pbn, 1))

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

test_that("generate_parametrisation_instance correctly generates instance", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  bn <- generate_parametrisation_instance(pbn, get_parametrisation_by_index(pbn, 2))

  expect_equivalent(bn$genes, c("G1", "G2", "G3"))

  expect_equivalent(bn$interactions$G1$input, c(3))
  expect_equivalent(bn$interactions$G1$func, c(0,1))

  expect_equivalent(bn$interactions$G2$input, c(1,2))
  expect_equivalent(bn$interactions$G2$func, c(1,0,1,1))

  expect_equivalent(bn$interactions$G3$input, c(1))
  expect_equivalent(bn$interactions$G3$func, c(0,1))

  expect_equivalent(bn$fixed, c(-1, -1, -1))
  expect_equivalent(class(bn), "BooleanNetwork")
})
