#' Computes the quantitative discrete bifurcation function of the
#' given baseline BN in the parametrisation space represented by the given RG.
#' If necessary, the motif extension is created before computation.
#' @param bn_baseline_network the baseline BN. Either a file path or a BooleanNetworkin-stance returned by load_BN.
#' @param rg_rewiring_space the RG representing the parametrisation space. Either a file path or an RG data frame.
#' @param output_genes the vector of output genes. Must be a subset of genes inbn. If not provided, defaults to BN genes
#' @param semantics  the BN semantics to be used:"async"or"sync"
#' @param attractor_similarity the attractor similarity measure to be used:"activity"or"overlap"
#' @param return_pbn the logical value indicating whether to also return the constructed PBN
#' @return The table of rewiring distance and attractor landscape similarity for all parametrisations.
#' @export
compute_discrete_bifurcation <- function(bn_baseline_network, rg_rewiring_space, output_genes = character(),
                                         semantics = "async", attractor_similarity = "activity", return_pbn = FALSE) {
  inputs <- process_inputs(bn_baseline_network, rg_rewiring_space, output_genes)

  intermediates <- create_extension(inputs$bn, inputs$rg, inputs$output_genes)

  parametrisations_count <- length(intermediates$pbn$function_index_combinations)

  message(paste("Starting computation for", parametrisations_count, "parametrisations."))

  if(parametrisations_count > 200000) {
    message(paste("Warning! This operation may take a significant time (hours)!"))
  } else if(parametrisations_count > 5000) {
    message(paste("Warning! This operation may take some time (minutes)!"))
  } else if(parametrisations_count > 1000) {
    message(paste("Warning! This operation may take some time (seconds)!"))
  }

  discrete_bifurcation <- compute_discrete_bifurcation_inner(intermediates$pbn, intermediates$bn_parametrisation,
                                                             output_genes = inputs$output_genes, semantics = semantics,
                                                             attractor_similarity = attractor_similarity)

  rewiring_distance <- sapply(discrete_bifurcation, function(row) row$rewiring_distance)
  attractor_landscape_similarity <- sapply(discrete_bifurcation, function(row) row$attractor_landscape_similarity)

  discrete_bifurcation_df <- data.frame(rewiring_distance, attractor_landscape_similarity)

  parameter_count <- sum(sapply(intermediates$pbn$gene_regulators_indexes, function(regulators_vector) 2^length(regulators_vector)))

  result <- list(discrete_bifurcation = discrete_bifurcation_df,
                 parameter_count = parameter_count)

  if(return_pbn) {
    result$pbn = intermediates$pbn
  }

  message("Operation finished successfully!")

  return(result)
}

#' The inner computation of QDBF.
compute_discrete_bifurcation_inner <- function(pbn, baseline_parametrisation, output_genes = character(),
                                               semantics = "async", attractor_similarity = "activity") {

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
                                           compute_db_single(index_combination_vector, gene_function_vectors,
                                                             baseline_parametrisation, baseline_attractor_landscape,
                                                             output_genes, output_genes_encoder,
                                                             genes, gene_regulators_indexes, instance_states, pre_to_post,
                                                             semantics, attractor_similarity),
                                         packages = c("BoolNet"), seq_threshold = 500,
                                         log_frequency = 500)

  return(discrete_bifurcation)
}

#' The QDBF copmutation for one given parametrisation
compute_db_single <- function(index_combination_vector, gene_function_vectors,
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

#' Computes the phenotype rewiring robustness over the data result-ing from a computation ofcompute_discrete_bifurcation.
#' @param discrete_bifurcation_result the quantitative discrete bifurcation returned by compute_discrete_bifurcation
#' @param p_elementary the elementary rewiring probability
#' @return The value of phenotype rewiring robustness
#' @export
compute_rewiring_robustness <- function(discrete_bifurcation_result, p_elementary = 0.1) {
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

#' Draws the bifurcation plot for the given data
#' @param discrete_bifurcation_result the quantitative discrete bifurcation returned by compute_discrete_bifurcation
#' @param main the title for the plot
#' @return The bifurcation plot
#' @export
bifurcation_plot <- function(discrete_bifurcation_result, main = "") {
  `%>%` <- magrittr::`%>%`

  distance_max <- discrete_bifurcation_result$parameter_count

  data <- discrete_bifurcation_result$discrete_bifurcation

  alpha_prior <- 16 / (nrow(discrete_bifurcation_result$discrete_bifurcation)^0.75)

  alpha <- max(1/256, min(1, alpha_prior))

  return(data %>%
    ggplot2::ggplot(ggplot2::aes(x=rewiring_distance, y=attractor_landscape_similarity, group = rewiring_distance,
               xmin = 0, xmax = distance_max, ymin = 0, ymax = 1)) +
    ggplot2::geom_boxplot(fill = "gray", outlier.size = 0.05) +
    ggplot2::geom_jitter(color="blue", alpha = alpha, width = 0.2, height = 0.01, size = 1.7, shape = 20) +
    ggplot2::ggtitle(main) +
    ggplot2::xlab("Rewiring distance") +
    ggplot2::ylab("Attractor landscape similarity"))
}

#' Get the attractor landscape of a given parametrisation in the form of attractors list
#' and attractor distribution matrix.
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

#' Process the inputs given to compute_discrete_bifurcation
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
