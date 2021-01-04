# O(ni^2 * 2^np)
#' Build the transition table in the form of an matrix:
#' row = index of the (all-genes) state (binary+1), col = (instance-gene)
build_transition_table <- function(instance, all_genes) {
  conversion_vectors <- lapply(instance$genes, (function(g) compute_conversion_vector(instance, all_genes, g)))

  input_rows <- all_binary_vectors(length(all_genes))

  transition_table <- matrix(0, nrow = length(input_rows), ncol = length(instance$genes))

  for (i in 1:length(input_rows)) {
    transition_table[i,] = compute_transition_table_row(instance, conversion_vectors, input_rows[[i]])
  }

  return(transition_table)
}

# O(ni * np)
#' Compute the conversion vector for input genes of a given gene.
#' The resulting vector contains indexes of input genes with regards to the all_genes vector.
compute_conversion_vector <- function(instance, all_genes, gene) {
  interaction <- instance$interactions[[gene]]
  result <- match(instance$genes[interaction$input], all_genes)

  return(result)
}

# O(ni^2)
#' Compute one row of transition table.
compute_transition_table_row <- function(instance, conversion_vectors, input_row) {
  cell_computer <- (function(gene) compute_transition_table_cell(instance, conversion_vectors, input_row, gene))
  transition_table_row <- sapply(instance$genes, cell_computer, USE.NAMES = FALSE)

  return(transition_table_row)
}

# O(ni)
#' Compute value of one cell of transition table, given input row and column gene.
compute_transition_table_cell <- function(instance, conversion_vectors, input_row, cell_gene) {
  gene_instance_index <- match(cell_gene, instance$genes)
  gene_regulators_vector <- conversion_vectors[[gene_instance_index]]
  gene_regulators_values <- input_row[gene_regulators_vector]
  gene_interaction_table_row_index <- binary_vector_to_one_based_index(gene_regulators_values)

  gene_interaction_table <- instance$interactions[[cell_gene]]$func
  gene_value_on_input <- gene_interaction_table[gene_interaction_table_row_index]

  return(gene_value_on_input)
}

#' Extract a BN instance from a PBN for a given index.
#' @param pbn The PBN to extract from
#' @param instance_index 1-based numeric index of a parametrisation in the PBN
#' @return the BN instance
#' @export
extract_instance <- function(pbn, instance_index) {
  parametrisation <- get_parametrisation_by_index(pbn, instance_index)

  instance <- generate_parametrisation_instance(pbn, parametrisation)

  return(instance)
}
