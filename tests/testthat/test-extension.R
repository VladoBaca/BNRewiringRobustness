context("Extension functionality of the package")

test_that("Extension is correctly created from both directions", {
  bn <- load_BN(system.file("examples", "motifs", "fb_minus.bn", package = "BNRewiringRobustness"))
  rg <- load_RG(system.file("examples", "rgs", "extension_test.rg", package = "BNRewiringRobustness"))
  output_genes <- c("A", "B")

  extension <- create_extension(bn, rg, output_genes)

  expect_equivalent(extension$pbn$genes, c("A", "B", "C", "E"))
  expect_equivalent(extension$pbn$gene_regulators_indexes, list(c(1,2,3), c(1), c(1,2,4), numeric()))

  expect_equivalent(extension$bn_parametrisation, list(c(1,1,0,0,1,1,0,0), c(0,1),  c(0,0,0,0,0,0,0,0), c(0)))

  bn <- generate_BN(extension$pbn$genes, extension$pbn$gene_regulators_indexes, extension$bn_parametrisation)

  expect_equivalent(bn$genes, c("A", "B", "C", "E"))

  expect_equivalent(bn$interactions$A$input, c(1,2,3))
  expect_equivalent(bn$interactions$A$func, c(1,1,0,0,1,1,0,0))

  expect_equivalent(bn$interactions$B$input, c(1))
  expect_equivalent(bn$interactions$B$func, c(0,1))

  expect_equivalent(bn$interactions$C$input, c(1,2,4))
  expect_equivalent(bn$interactions$C$func, c(0,0,0,0,0,0,0,0))

  expect_equivalent(bn$interactions$E$input, c(0))
  expect_equivalent(bn$interactions$E$func, c(0))
})

test_that("create_extension correctly validates genes and edges", {
  rg_A <- load_RG(system.file("examples", "rgs", "general_1.rg", package = "BNRewiringRobustness"))
  rg_ABCDE_sparse <- load_RG(system.file("examples", "rgs", "extension_test.rg", package = "BNRewiringRobustness"))

  bn_A <- load_BN(system.file("examples", "motifs", "auto_plus.bn", package = "BNRewiringRobustness"))
  bn_AB_dense <- load_BN(system.file("examples", "motifs", "crm_and.bn", package = "BNRewiringRobustness"))

  expect_error({create_extension(bn_A, rg_A, c("A", "B"))}, regexp = "Output genes must be a subset of BN genes!")
  expect_error({create_extension(bn_AB_dense, rg_A, c("A"))}, regexp = "BN genes must be a subset of RG genes!")
  expect_error({create_extension(bn_AB_dense, rg_ABCDE_sparse, c("A", "B"))}, regexp = "BN edges must be a subset of RG edges!")
})
