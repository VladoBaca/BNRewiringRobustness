context("IO functionality of the package")

test_that("load_BN correctly loads the Boolean networks", {
  fb_minus_net <- load_BN(system.file("examples", "motifs", "fb_minus.bn", package = "BNRewiringRobustness"))

  expect_equivalent(fb_minus_net$genes, c("A", "B"))
  expect_equivalent(fb_minus_net$fixed, c(-1, -1))
  expect_equivalent(fb_minus_net$interactions$A$input, c(2))
  expect_equivalent(fb_minus_net$interactions$A$func, c(1, 0))
  expect_equivalent(fb_minus_net$interactions$B$input, c(1))
  expect_equivalent(fb_minus_net$interactions$B$func, c(0, 1))

  ffl_c_1 <- load_BN(system.file("examples", "motifs", "ffl_c1.bn", package = "BNRewiringRobustness"))

  expect_equivalent(ffl_c_1$genes, c("A", "B", "C"))
  expect_equivalent(ffl_c_1$fixed, c(-1, -1, -1))
  expect_equivalent(ffl_c_1$interactions$A$input, c(1))
  expect_equivalent(ffl_c_1$interactions$A$func, c(1, 1))
  expect_equivalent(ffl_c_1$interactions$B$input, c(1))
  expect_equivalent(ffl_c_1$interactions$B$func, c(0, 1))
  expect_equivalent(ffl_c_1$interactions$C$input, c(1, 2))
  expect_equivalent(ffl_c_1$interactions$C$func, c(0, 1, 1, 1))

  crm_assymetric <- load_BN(system.file("examples", "motifs", "crm_assymetric.bn", package = "BNRewiringRobustness"))

  expect_equivalent(crm_assymetric$genes, c("A", "B"))
  expect_equivalent(crm_assymetric$fixed, c(-1, -1))
  expect_equivalent(crm_assymetric$interactions$A$input, c(1, 2))
  expect_equivalent(crm_assymetric$interactions$A$func, c(1, 0, 1, 1))
  expect_equivalent(crm_assymetric$interactions$B$input, c(1, 2))
  expect_equivalent(crm_assymetric$interactions$B$func, c(0, 1, 0, 0))
})


test_that("load_PBN correctly loads the parametrised Boolean network", {
  test_net <- load_PBN(system.file("examples", "rgs", "all_reg_types.rg", package = "BNRewiringRobustness"))

  expect_equivalent(test_net$genes, c("A", "B", "C", "D"))

  expected_function_vectors <- list(
    list(c(1,0,0,0),c(1,1,0,0),c(1,1,1,0)),
    list(c(0,0,0,1),c(0,0,1,0),c(0,0,1,1),c(0,1,1,1),c(1,0,1,1)),
    list(c(0,0,0,1),c(0,0,1,1),c(0,1,0,0),c(0,1,1,1),c(1,1,0,0),c(1,1,0,1)),
    list(c(0,1))
  )

  expect_equivalent(test_net$gene_function_vectors, expected_function_vectors)

  expect_equivalent(test_net$gene_regulators_indexes, list(c(2,3),c(1,3),c(1,2),c(4)))

  combs <- expand.grid(1:3, 1:5, 1:6, 1:1)
  m <- apply(combs, 2, (function(x) x))
  expected_combinations <- (unname(tapply(m,rep(1:nrow(m),ncol(m)),function(i)i)))

  expect_equivalent(test_net$function_index_combinations,expected_combinations)
})

test_that("load_RG correctly loads the regulatory graph", {
  test_graph <- load_RG(system.file("examples", "rgs", "test_middle.rg", package = "BNRewiringRobustness"))

  n1 <- c("A", "B","C","A","B","C","A","B","C")
  edge <- c("-?","-?","-??" ,"-?","->?" ,"-?","-|?", "-?", "->")
  n2 <- c("A", "A","A","B","B","B","C","C","C")
  expected <- data.frame(n1, edge, n2)

  expect_equivalent(test_graph, expected)
})
