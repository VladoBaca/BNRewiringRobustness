#' @export
compute_discrete_bifurcation <- function(bn_baseline_network, rg_rewiring_space, output_genes = character(),
                                         semantics = "async", attractor_similarity = "activity", return_pbn = FALSE, verbose = FALSE) {
  inputs <- process_inputs(bn_baseline_network, rg_rewiring_space, output_genes)

  intermediates <- create_extension(inputs$bn, inputs$rg, inputs$output_genes)

  discrete_bifurcation <- compute_discrete_bifurcation_inner(intermediates$pbn, intermediates$bn_parametrisation,
                                                             output_genes = inputs$output_genes, semantics = semantics,
                                                             attractor_similarity = attractor_similarity, verbose)

  rewiring_distance <- sapply(discrete_bifurcation, function(row) row$rewiring_distance)
  attractor_landscape_similarity <- sapply(discrete_bifurcation, function(row) row$attractor_landscape_similarity)

  discrete_bifurcation_df <- data.frame(rewiring_distance, attractor_landscape_similarity)

  parameter_count <- sum(sapply(intermediates$pbn$gene_regulators_indexes, function(regulators_vector) 2^length(regulators_vector)))

  result <- list(discrete_bifurcation = discrete_bifurcation_df,
                 parameter_count = parameter_count)

  if(return_pbn) {
    result$pbn = intermediates$pbn
  }

  return(result)
}

compute_discrete_bifurcation_inner <- function(pbn, baseline_parametrisation, output_genes = character(),
                                               semantics = "async", attractor_similarity = "activity", verbose = FALSE) {

  genes <- pbn$genes
  instance_states <- all_binary_vectors(length(genes))
  pre_to_post <- pre_to_post_vector(instance_states)
  output_genes_encoder <- match(output_genes, genes)
  baseline_attractor_landscape <- get_attractor_landscape(baseline_parametrisation, genes, pbn$gene_regulators_indexes,
                                                          instance_states, pre_to_post, semantics)

  gene_regulators_indexes <- pbn$gene_regulators_indexes

  gene_function_vectors <- pbn$gene_function_vectors

  discrete_bifurcation <- parallel_apply(pbn$function_index_combinations,
                                         function(index_combination_vector)
                                           compute_db_single(index_combination_vector, gene_function_vectors,  #parametrisation,
                                                             baseline_parametrisation, baseline_attractor_landscape,
                                                             output_genes, output_genes_encoder,
                                                             genes, gene_regulators_indexes, instance_states, pre_to_post,
                                                             semantics, attractor_similarity),
                                           packages = c("BoolNet"), seq_threshold = 500,
                                           message_frequency = 500, verbose = verbose)

  return(discrete_bifurcation)
}

compute_db_single <- function(index_combination_vector, gene_function_vectors, #parametrisation,
                              baseline_parametrisation, baseline_attractor_landscape,
                              output_genes, output_genes_encoder,
                              genes, gene_regulators_indexes, instance_states, pre_to_post,
                              semantics, attractor_similarity) {
  parametrisation <- get_parametrisation_by_index_combination_vector(gene_function_vectors, index_combination_vector)

  rewiring_distance <- compute_rewiring_distance(baseline_parametrisation, parametrisation)

  attractor_landscape <- get_attractor_landscape(parametrisation, genes, gene_regulators_indexes,
                                                 instance_states, pre_to_post, semantics)

  attractor_landscape_similarity <- compute_attractor_landscape_similarity(baseline_attractor_landscape, attractor_landscape,
                                                                           genes, output_genes, output_genes_encoder, instance_states, attractor_similarity)

  return(list(rewiring_distance = rewiring_distance, attractor_landscape_similarity = attractor_landscape_similarity))
}

#TODO the default prob
#' @export
compute_rewiring_robustness <- function(discrete_bifurcation_result, p_elementary = 0.001) {
  parametrisation_count <- nrow(discrete_bifurcation_result$discrete_bifurcation)

  probabilities_sum <- 0
  robustness <- 0

  for (parametrisation_i in 1:parametrisation_count) {
    distance <- discrete_bifurcation_result$discrete_bifurcation$rewiring_distance[parametrisation_i]
    similarity <- discrete_bifurcation_result$discrete_bifurcation$attractor_landscape_similarity[parametrisation_i]

    probability_prior <- p_elementary^distance * (1 - p_elementary)^(discrete_bifurcation_result$parameter_count - distance)

    probabilities_sum <- probabilities_sum + probability_prior
    robustness <- robustness + probability_prior * similarity
  }

  robustness <- robustness / probabilities_sum
  return(robustness)
}

#TODO experiment with this?
#' @export
bifurcation_plot <- function(discrete_bifurcation_result) {
  distance_max <- discrete_bifurcation_result$parameter_count

  distance <- discrete_bifurcation_result$discrete_bifurcation$rewiring_distance
  attractor_landscape_similarity <- discrete_bifurcation_result$discrete_bifurcation$attractor_landscape_similarity

  plot(distance, attractor_landscape_similarity,
       xlim=c(0, distance_max), ylim=c(0, 1),
       xlab = "Rewiring distance", ylab = "Attractor landscape similarity",
       col = rgb(red = 0, green = 0, blue = 1, alpha = 0.2))

  #abline(lm(attractor_landscape_similarity ~ distance))
}

# Get the attractor landscape of a given parametrisation in the form of attractors list
# and attractor distribution matrix.
get_attractor_landscape <- function(parametrisation, genes, gene_regulators_indexes, instance_states, pre_to_post, semantics = "async")
{
  bn <- generate_BN(genes, gene_regulators_indexes, parametrisation)

  if (semantics == "sync" | semantics == "synchronous") {
    instance_attractors_info <- get_attractors_data_sync(bn, pre_to_post)
  } else {
    instance_attractors_info <- get_attractors_data_async(bn, instance_states, pre_to_post)
  }

  return(instance_attractors_info)
}

process_inputs <- function(bn_baseline_network, rg_rewiring_space, output_genes) {
  if (is_simple_string(bn_baseline_network)) {
    bn <- load_BN(bn_baseline_network)
  } else {
    bn <- bn_baseline_network
  }

  if (is_simple_string(rg_rewiring_space)) {
    rg <- load_RG(rg_rewiring_space)
  } else {
    rg <- rg_rewiring_space
  }

  if(length(output_genes) == 0) {
    output_genes <- bn$genes
  }

  return(list(bn = bn, rg = rg, output_genes = output_genes))
}
