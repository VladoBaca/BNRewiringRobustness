compute_rewiring_distance <- function(parametrisation_1, parametrisation_2) {
  distances_per_gene <- sapply(1:length(parametrisation_1),
                                  function(gene_i) sum(abs(parametrisation_1[[gene_i]] - parametrisation_2[[gene_i]])))
  distance <- sum(distances_per_gene)

  return(distance)
}

compute_attractor_landscape_similarity <- function(attractor_landscape_1, attractor_landscape_2,
                                                   genes, output_genes, output_genes_encoder, instance_states,
                                                   attractor_similarity = "activity") {

  attractors_similarity_matrix <- compute_attractors_similarity_matrix(attractor_landscape_1$attractors,
                                                                       attractor_landscape_2$attractors,
                                                                       genes, output_genes, output_genes_encoder,
                                                                       instance_states, attractor_similarity)

  lp_objective_coefficients <- as.vector(t(attractors_similarity_matrix))
  lp_constraints_coefficients_matrix <- compute_constraints_coefficients_matrix(length(attractor_landscape_1$attractors),
                                                                                length(attractor_landscape_2$attractors))
  lp_directions <- rep("==", length(attractor_landscape_1$attractors) + length(attractor_landscape_2$attractors))

  attractor_distribution_similarity_vector <- sapply(1:length(instance_states), function(state_i)
    compute_attractor_distribution_similarity(attractor_landscape_1$distribution_matrix[state_i,],
                                              attractor_landscape_2$distribution_matrix[state_i,],
                                              lp_objective_coefficients, lp_constraints_coefficients_matrix, lp_directions))

  attractor_landscape_similarity <- mean(attractor_distribution_similarity_vector)
  return(attractor_landscape_similarity)
}

compute_attractors_similarity_matrix <- function(attractors_1, attractors_2, genes, output_genes, output_genes_encoder,
                                                 instance_states, attractor_similarity) {
  if (attractor_similarity == "overlap") {
    return(compute_attractors_similarity_matrix_overlap(attractors_1, attractors_2, output_genes,
                                                        output_genes_encoder, instance_states))
  } else {
    return(compute_attractors_similarity_matrix_activity(attractors_1, attractors_2, genes,
                                                         output_genes, output_genes_encoder, instance_states))
  }
}

compute_attractors_similarity_matrix_activity <- function(attractors_1, attractors_2, genes,
                                                          output_genes, output_genes_encoder, instance_states) {
  activity_ratios_matrix_1 <- compute_activity_ratios_matrix(attractors_1, genes, output_genes,
                                                             output_genes_encoder, instance_states)
  activity_ratios_matrix_2 <- compute_activity_ratios_matrix(attractors_2, genes, output_genes,
                                                             output_genes_encoder, instance_states)

  attractors_similarity_matrix <- matrix(0, nrow = length(attractors_1), ncol = length(attractors_2))

  for (i in 1:length(attractors_1)) {
    for (j in 1:length(attractors_2)) {
      attractors_similarity_matrix[i,j] <- sum(1 - abs(activity_ratios_matrix_1[,i] - activity_ratios_matrix_2[,j])) / length(output_genes)
    }
  }

  return(attractors_similarity_matrix)
}

compute_activity_ratios_matrix <- function(attractors, genes, output_genes, output_genes_encoder, instance_states) {
  attractors_activity_ratios_matrix <- matrix(0, nrow = length(output_genes), ncol = length(attractors))

  output_states <- lapply(instance_states, function(state_vector) state_vector[output_genes_encoder])

  for (attractor_i in 1:length(attractors)) {
    attractor_activity_ratios <- rep(0, length(output_genes))

    for (state_i0 in attractors[[attractor_i]]) {
      attractor_activity_ratios <- attractor_activity_ratios + output_states[[state_i0 + 1]]
    }

    attractor_activity_ratios <- attractor_activity_ratios / length(attractors[[attractor_i]])

    attractors_activity_ratios_matrix[,attractor_i] <- attractor_activity_ratios
  }
  return(attractors_activity_ratios_matrix)
}

compute_attractors_similarity_matrix_overlap <- function(attractors_1, attractors_2, output_genes,
                                                         output_genes_encoder, instance_states) {
  attractors_output_states_ratios_matrix_1 <- compute_attractors_output_states_ratios_matrix(attractors_1, output_genes,
                                                                                             output_genes_encoder,
                                                                                             instance_states)
  attractors_output_states_ratios_matrix_2 <- compute_attractors_output_states_ratios_matrix(attractors_2, output_genes,
                                                                                             output_genes_encoder,
                                                                                             instance_states)

  attractors_similarity_matrix <- matrix(0, nrow = length(attractors_1), ncol = length(attractors_2))

  output_states_count <- 2^length(output_genes)

  for (attractor_1_i in 1:length(attractors_1)) {
    for (attractor_2_i in 1:length(attractors_2)) {
      for (output_state_i in 1:output_states_count) {
        overlap <- min(attractors_output_states_ratios_matrix_1[output_state_i, attractor_1_i],
                       attractors_output_states_ratios_matrix_2[output_state_i, attractor_2_i])

        attractors_similarity_matrix[attractor_1_i, attractor_2_i] <-
          attractors_similarity_matrix[attractor_1_i, attractor_2_i] + overlap
      }
    }
  }

  return(attractors_similarity_matrix)
}

compute_attractors_output_states_ratios_matrix <- function(attractors, output_genes,
                                                           output_genes_encoder, instance_states) {
  output_states_count <- 2^length(output_genes)

  attractors_output_states_ratios_matrix <- matrix(0, nrow = output_states_count, ncol = length(attractors))

  for (attractor_i in 1:length(attractors)) {
    state_unit_weight <- 1 / length(attractors[[attractor_i]])
    for (state_i0 in attractors[[attractor_i]]) {
      output_state_index <- binary_vector_to_one_based_index(instance_states[[state_i0 + 1]][output_genes_encoder])

      attractors_output_states_ratios_matrix[output_state_index,attractor_i] <-
        attractors_output_states_ratios_matrix[output_state_index,attractor_i] + state_unit_weight
    }
  }

  return(attractors_output_states_ratios_matrix)
}

compute_constraints_coefficients_matrix <- function(attractor_1_count, attractor_2_count) {
  constraints_coefficients_matrix <- matrix(0, nrow = attractor_1_count + attractor_2_count, ncol = attractor_1_count * attractor_2_count)
  ones_1 <- rep(1, attractor_2_count)
  ones_2 <- rep(1, attractor_1_count)

  for (i in 1:attractor_1_count) {
    cols_vector <- ((i-1)*attractor_2_count + 1) : (i*attractor_2_count)
    constraints_coefficients_matrix[i, cols_vector] <- ones_1
  }

  for (j in 1:attractor_2_count) {
    cols_vector <- ((0 : (attractor_1_count - 1)) * attractor_2_count) + j
    constraints_coefficients_matrix[attractor_1_count + j, cols_vector] <- ones_2
  }

  return(constraints_coefficients_matrix)
}

compute_attractor_distribution_similarity <- function(attractor_distribution_1, attractor_distribution_2,
                                                      lp_objective_coefficients, lp_constraints_coefficients_matrix, lp_directions) {
  lp_constraint_sums <- c(attractor_distribution_1, attractor_distribution_2)

  solution <- lpSolve::lp(direction = "max", lp_objective_coefficients,
                          lp_constraints_coefficients_matrix, lp_directions, lp_constraint_sums)

  return(solution$objval)
}
