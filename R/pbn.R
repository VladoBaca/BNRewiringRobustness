#' Creates a PBN from given RG.
#' @rg The regulatory graph
#' @return The PBN.
#' @export
create_pbn_from_rg <- function(rg) {
  validate_unique_edges(rg)

  genes <- get_rg_genes(rg)

  gene_function_vectors_full <- parallel_apply(genes, function(gene) {enumerate_gene_function_vectors(rg, gene)}, seq_threshold = 5)

  gene_regulators_indexes <- lapply(gene_function_vectors_full, function(r) match(r$regulators, genes))

  gene_function_vectors <- lapply(gene_function_vectors_full, function(r) r$function_vectors)

  function_index_combinations <- as.matrix(expand.grid(lapply(gene_function_vectors, function(p) 1:length(p))))

  function_index_combinations <- unname(tapply(function_index_combinations,
                                               rep(1:nrow(function_index_combinations),
                                                   ncol(function_index_combinations)),function(i)i))

  return(list(genes = genes, gene_regulators_indexes = gene_regulators_indexes,
              gene_function_vectors = gene_function_vectors, function_index_combinations = function_index_combinations))
}

#' Validate that there are no repeating edges in the RG
validate_unique_edges <- function(rg) {
  all_unique <- all(sapply(1:nrow(rg), function(i) {length(which(rg$n1 == rg$n1[i] & rg$n2 == rg$n2[i])) == 1 } ))

  if(!all_unique) {
    stop("Not all the edges are unique in the regulatory graph!")
  }
}

#' Get all the genes of the given RG.
get_rg_genes <- function(rg) {
  return(union(levels(as.factor(rg$n1)), levels(as.factor(rg$n2))))
}

# O(2^(2^n) * 2^n * n)
#' Enumerates all valid function vectors for given gene
enumerate_gene_function_vectors <- function(rg, gene) {
  regulation_rows <- rg[which(rg$n2 == gene),]
  regulation_count <- nrow(regulation_rows)

  regulators <- sapply(regulation_rows$n1, as.character)
  regulations <- lapply(regulation_rows$edge, as.character)

  row_count <- 2^regulation_count

  all_possible_function_vectors <- all_binary_vectors(row_count)

  input_rows <- all_binary_vectors(regulation_count)

  function_vector_validator <- (function(function_vector) is_function_vector_valid(function_vector, regulation_count, regulations, input_rows))

  valid_function_vectors <- all_possible_function_vectors[sapply(all_possible_function_vectors, function_vector_validator)]

  return(list(regulators = regulators, function_vectors = valid_function_vectors))
}

# O(n 2^n)
#' Validates that given function vector satisfies the constraints given by the regulatory graph.
is_function_vector_valid <- function(function_vector, regulation_count, regulations, input_rows) {
  if(regulation_count == 0)
    return(TRUE)

  regulation_observable_vector <- rep(FALSE, regulation_count)

  for (input_row_index in 0:(length(input_rows) - 1)) {
    input_row <- input_rows[[input_row_index + 1]]

    for (input_variable_index in 0:(regulation_count - 1)) {
      input_variable <- input_row[input_variable_index + 1]

      if (input_variable == 0) {
        flipped_index <- input_row_index + 2^(regulation_count - input_variable_index - 1)
        output_value <- function_vector[input_row_index + 1]
        flipped_value <- function_vector[flipped_index + 1]

        if (output_value != flipped_value) {
          regulation_observable_vector[input_variable_index + 1] <- TRUE
        }

        regulation <- regulations[[input_variable_index + 1]]
        if ((regulation == "->" || regulation == "->?") && flipped_value < output_value)
          return(FALSE)
        if ((regulation == "-|" || regulation == "-|?") && flipped_value > output_value)
          return(FALSE)
      }
    }
  }

  for (regulation_index in 0:(regulation_count - 1)) {
    regulation <- regulations[[regulation_index + 1]]
    if ((regulation == "->" || regulation == "-|" || regulation == "-?") && !regulation_observable_vector[regulation_index + 1]) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Generates a BooleanNetwork instance from the given data
generate_BN <- function(genes, gene_regulators_indexes, gene_function_vectors) {
  fixed <- rep(-1, length(genes))
  names(fixed) <- genes

  interactions <- lapply(1:length(genes), function(gene_i) generate_interaction(gene_i, gene_regulators_indexes, gene_function_vectors))
  names(interactions) <- genes

  network <- list(interactions = interactions, genes = genes, fixed = fixed)
  class(network) <- "BooleanNetwork"
  return(network)
}

#' Generates a single interaction for genereated BN
generate_interaction <- function(gene_i, gene_regulators_indexes, gene_function_vectors) {
  if (length(gene_regulators_indexes[[gene_i]]) == 0) {
    input <- 0
  } else {
    input <- gene_regulators_indexes[[gene_i]]
  }
  return(list(input = input, func = gene_function_vectors[[gene_i]], expression = ""))
}

#' Create a BN from a given parametrisation
generate_parametrisation_instance <- function(pbn, parametrisation) {
  instance <- generate_BN(pbn$genes, pbn$gene_regulators_indexes, parametrisation)

  return(instance)
}

#' Get a parametrisation from PBN by its numeric 1-based index.
get_parametrisation_by_index <- function(pbn, parametrisation_index) {
  index_combination_vector <- pbn$function_index_combinations[[parametrisation_index]]  # as.numeric(pbn$function_index_combinations[parametrisation_index,])

  return(get_parametrisation_by_index_combination_vector(pbn$gene_function_vectors, index_combination_vector))
}

#' Get a parametrisation from PBN by vector of indexes of update functions of respective genes.
get_parametrisation_by_index_combination_vector <- function(gene_function_vectors, index_combination_vector) {
  parametrisation <- lapply(1:length(index_combination_vector),
                            function(i) gene_function_vectors[[i]][[index_combination_vector[i]]])
  return(parametrisation)
}
