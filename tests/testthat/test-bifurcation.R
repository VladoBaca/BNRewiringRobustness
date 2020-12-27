context("Functionality of computing the quantitative bifurcation and robustness")

test_that("compute_discrete_bifurcation computes correct quantitative bifurcation function", {
  rg <- load_RG(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  bn <- load_BN(system.file("examples", "motifs", "test_triple_and.bn", package = "BNRewiringRobustness"))

  output_genes = c("G1", "G2")

  rewiring_distance <- c(0, 2)
  attractor_landscape_similarity <- c(1, 0.75)

  expected_result <- list(discrete_bifurcation = data.frame(rewiring_distance, attractor_landscape_similarity),
                          parameter_count = 8)

  expected_result_pbn <- list(discrete_bifurcation = data.frame(rewiring_distance, attractor_landscape_similarity),
                          parameter_count = 8, pbn = create_pbn_from_rg(rg))


  actual_result_no_pbn <- compute_discrete_bifurcation(bn, rg, output_genes = output_genes,
                                           semantics = "async", attractor_similarity = "activity", return_pbn = FALSE)

  actual_filenames_pbn <- compute_discrete_bifurcation(system.file("examples", "motifs", "test_triple_and.bn", package = "BNRewiringRobustness"),
                                                       load_RG(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness")),
                                                       output_genes = output_genes,
                                                       semantics = "async", attractor_similarity = "activity", return_pbn = TRUE)

  expect_equivalent(actual_result_no_pbn, expected_result)
  expect_equivalent(actual_filenames_pbn, expected_result_pbn)
})

test_that("compute_rewiring_robustness computes correct rewiring robustness", {
  rg <- load_RG(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))
  bn <- load_BN(system.file("examples", "motifs", "test_triple_and.bn", package = "BNRewiringRobustness"))

  output_genes = c("G1", "G2")

  qdb <- compute_discrete_bifurcation(bn, rg, output_genes = output_genes,
                                      semantics = "async", attractor_similarity = "activity", return_pbn = FALSE)

  actual_robustness <- compute_rewiring_robustness(qdb, p_elementary = 0.1)

  prior_prob_sum <- 0.43578162

  psi_1 <- 0.43046721 / prior_prob_sum
  psi_2 <- 0.00531441 / prior_prob_sum

  expected_robustness <- psi_1 * 1 + psi_2 * 0.75

  expect_equivalent(actual_robustness, expected_robustness)
})
