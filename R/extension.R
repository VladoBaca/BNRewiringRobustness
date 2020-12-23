#' Transform the given motif BN and possibly larger PBN to an intermediate gene set.
create_extension <- function(bn_motif, rg_large, output_genes) {
  validate_genes(bn_motif, rg_large, output_genes)
  validate_edges(bn_motif, rg_large)

  intermediate_genes <- compute_intermediate_genes(bn_motif$genes, rg_large)

  rg_intermediate <- create_intermediate_rg(rg_large, intermediate_genes)

  pbn_intermediate <- create_pbn_from_rg(rg_intermediate)

  bn_intermediate_parametrisation <- create_intermediate_motif_parametrisation(bn_motif, pbn_intermediate)

  return(list(pbn = pbn_intermediate, bn_parametrisation = bn_intermediate_parametrisation))
}

#' Compute the set of transitive regulators of motif genes
compute_intermediate_genes <- function(motif_genes, rg_large)
{
  gene_set <- motif_genes
  changed <- TRUE

  while (changed) {
    regulators <- unique(rg_large[rg_large$n2 %in% gene_set, "n1"])
    gene_union <- union(gene_set, regulators)
    changed <- length(gene_union) > length(gene_set)
    gene_set <- gene_union
  }

  return(gene_set)
}

create_intermediate_rg <- function(rg, intermediate_genes) {
  intermediate_rg <- rg[rg$n2 %in% intermediate_genes,]
  return(intermediate_rg)
}

create_intermediate_motif_parametrisation <- function(bn_motif, pbn_intermediate) {
  motif_transition_table <- build_transition_table(bn_motif, pbn_intermediate$genes)

  parametrisation <- lapply(1:length(pbn_intermediate$genes),
                                  function(gene_i) get_function_vector(gene_i, bn_motif, pbn_intermediate, motif_transition_table))
  return(parametrisation)
}

#TODO this could be done without building transition table
get_function_vector <- function(gene_i, bn_motif, pbn_intermediate, motif_transition_table) {
  gene <- pbn_intermediate$genes[gene_i]
  regulators_count <- length(pbn_intermediate$gene_regulators_indexes[[gene_i]])
  function_vector_length <- 2^regulators_count

  function_vector <- rep(0,function_vector_length)

  if (gene %in% bn_motif$genes) {
    regulators_values_vectors <- all_binary_vectors(regulators_count)

    table_col_index <- match(gene, bn_motif$genes)

    for (function_vector_index in 1:function_vector_length) {
      function_vector_input_row <- regulators_values_vectors[[function_vector_index]]

      table_input_row <- rep(0,length(pbn_intermediate$genes))

      table_input_row[pbn_intermediate$gene_regulators_indexes[[gene_i]]] <- function_vector_input_row
      function_vector[function_vector_index] <- motif_transition_table[binary_vector_to_one_based_index(table_input_row), table_col_index]
    }
  }

  return(function_vector)
}

validate_genes <- function(bn_motif, rg_large, output_genes)
{
  if (!is_subset(output_genes, bn_motif$genes))
    stop("Output genes must be a subset of BN genes!")

  rg_large_genes <- get_rg_genes(rg_large)

  if (!is_subset(bn_motif$genes, rg_large_genes))
    stop("BN genes must be a subset of RG genes!")
}

validate_edges <- function(bn_motif, rg_large)
{
  for (gene_i in 1:length(bn_motif$genes)) {
    if(!is_function_vector_constant(bn_motif$interactions[[gene_i]]$func)) {
      bn_node_regulators <- bn_motif$genes[bn_motif$interactions[[gene_i]]$input]
      rg_node_regulators <- unique(rg_large[rg_large$n2 == bn_motif$genes[gene_i], "n1"])

      if (!is_subset(bn_node_regulators, rg_node_regulators))
        stop("BN edges must be a subset of RG edges!")
    }
  }
}

is_subset <- function(subset, superset) {
  genes_union <- union(superset, subset)

  return(length(genes_union) == length(superset))
}

is_function_vector_constant <- function(function_output_vector) {
  return(all(function_output_vector == 0) | all(function_output_vector == 1))
}


