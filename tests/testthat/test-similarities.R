context("Functionality of computing the different similarity measures")

test_that("compute_rewiring_distance computes correct distance", {
  pbn <- load_PBN(system.file("examples", "rgs", "test_triple.rg", package = "BNRewiringRobustness"))

  parametrisation_and <- get_parametrisation_by_index(pbn, 1)
  parametrisation_or <- get_parametrisation_by_index(pbn, 2)

  expect_equivalent(compute_rewiring_distance(parametrisation_and, parametrisation_and), 0)
  expect_equivalent(compute_rewiring_distance(parametrisation_and, parametrisation_or), 2)
})

test_that("compute_attractors_similarity_matrix_overlap computes correct similiarity", {
  attractors_1 <- list(c(7),c(0,2))
  attractors_2 <- list(c(0),c(5,7))

  output_genes <- c("G1", "G2")
  output_genes_encoder <- c(1,2)
  instance_states <- all_binary_vectors(3)

  sim_matrix <- compute_attractors_similarity_matrix_overlap(attractors_1, attractors_2, output_genes,
                                                             output_genes_encoder, instance_states)

  expected_matrix <- matrix(c(0, 0.5, 0.5, 0), nrow = 2, ncol = 2, byrow = T)

  expect_equivalent(sim_matrix, expected_matrix)
})

test_that("compute_activity_ratios_matrix computes correct activity ratios", {
  attractors <- list(c(7),c(0,2))

  genes <- c("G1", "G2", "G3")
  output_genes <- c("G1", "G2")
  output_genes_encoder <- c(1,2)
  instance_states <- all_binary_vectors(3)

  act_matrix <- compute_activity_ratios_matrix(attractors, genes, output_genes, output_genes_encoder, instance_states)

  expected_matrix <- matrix(c(1, 0, 1, 0.5), nrow = 2, ncol = 2, byrow = T)

  expect_equivalent(act_matrix, expected_matrix)
})

test_that("compute_attractors_similarity_matrix_activity computes correct similiarity", {
  attractors_1 <- list(c(7),c(0,2))
  attractors_2 <- list(c(0),c(5,7))

  genes <- c("G1", "G2", "G3")
  output_genes <- c("G1", "G2")
  output_genes_encoder <- c(1,2)
  instance_states <- all_binary_vectors(3)

  sim_matrix <- compute_attractors_similarity_matrix_activity(attractors_1, attractors_2, genes,
                                                              output_genes, output_genes_encoder, instance_states)

  expected_matrix <- matrix(c(0, 0.75, 0.75, 0.5), nrow = 2, ncol = 2, byrow = T)

  expect_equivalent(sim_matrix, expected_matrix)
})

test_that("compute_attractor_distribution_similarity computes correct similiarity", {
  attractors_1 <- list(c(7),c(0,2))
  attractors_2 <- list(c(0),c(5,7))

  attractors_similarity_matrix <- matrix(c(0, 0.75, 0.75, 0.5), nrow = 2, ncol = 2, byrow = T)

  genes <- c("G1", "G2", "G3")
  output_genes <- c("G1", "G2")
  output_genes_encoder <- c(1,2)
  instance_states <- all_binary_vectors(3)

  lp_objective_coefficients <- as.vector(t(attractors_similarity_matrix))
  lp_constraints_coefficients_matrix <- compute_constraints_coefficients_matrix(2, 2)
  lp_directions <- rep("==", 4)

  sim_1 <- compute_attractor_distribution_similarity(c(0, 1), c(1, 0),
                                                     lp_objective_coefficients, lp_constraints_coefficients_matrix, lp_directions)
  sim_2 <- compute_attractor_distribution_similarity(c(1, 0), c(0, 1),
                                                     lp_objective_coefficients, lp_constraints_coefficients_matrix, lp_directions)
  sim_3 <- compute_attractor_distribution_similarity(c(0.5, 0.5), c(0.5, 0.5),
                                                     lp_objective_coefficients, lp_constraints_coefficients_matrix, lp_directions)

  expect_equivalent(sim_1, 0.75)
  expect_equivalent(sim_2, 0.75)
  expect_equivalent(sim_3, 0.75)
})

test_that("compute_attractor_landscape_similarity computes correct similiarity", {
  attractors_1 <- list(c(7),c(0,2))
  matrix_1 <- matrix(c(0, 1,
                       0.5,0.5,
                       0, 1,
                       0.5,0.5,
                       0.5,0.5,
                       1, 0,
                       0.5,0.5,
                       1, 0), nrow = 8, ncol = 2, byrow = TRUE)

  attractors_2 <- list(c(0),c(5,7))
  matrix_2 <- matrix(c(1, 0,
                       0.5,0.5,
                       1, 0,
                       0.5,0.5,
                       0.5,0.5,
                       0, 1,
                       0.5,0.5,
                       0, 1), nrow = 8, ncol = 2, byrow = TRUE)

  attractor_landscape_1 <- list(attractors = attractors_1, distribution_matrix = matrix_1)
  attractor_landscape_2 <- list(attractors = attractors_2, distribution_matrix = matrix_2)

  genes <- c("G1", "G2", "G3")
  output_genes <- c("G1", "G2")
  output_genes_encoder <- c(1,2)
  instance_states <- all_binary_vectors(3)

  sim_activity <- compute_attractor_landscape_similarity(attractor_landscape_1, attractor_landscape_2,
                                                  genes, output_genes, output_genes_encoder, instance_states,
                                                  attractor_similarity = "activity")

  sim_overlap <- compute_attractor_landscape_similarity(attractor_landscape_1, attractor_landscape_2,
                                                        genes, output_genes, output_genes_encoder, instance_states,
                                                        attractor_similarity = "overlap")

  expect_equivalent(sim_activity, 0.75)
  expect_equivalent(sim_overlap, 0.5)
})

