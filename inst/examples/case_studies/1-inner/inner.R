library("BNRewiringRobustness")

# Prepare the data
rg <- load_RG(system.file("examples", "rgs", "general_2.rg", package = "BNRewiringRobustness"))

motif_names <- c("fb_plus", "fb_minus", "fb_double_minus", "crm_and", "crm_assymetric", "crm_or")
motif_nice_names <- c("DPFB", "NFB", "DNFB", "CRM (and)", "CRM (assymetric)", "CRM (or)")

motifs <- lapply(motif_names, function(motif_name)
  load_BN(system.file("examples", "motifs", paste0(motif_name, ".bn"), package = "BNRewiringRobustness")))

# Compute the discrete bifurcation
motif_qdbs <- lapply(motifs, function(motif) compute_discrete_bifurcation(motif, rg, return_pbn = TRUE))
pbn <- motif_qdbs[[1]]$pbn

read_result_from_file <- function(i) {
  discrete_bifurcation <- read.csv(paste0(motif_names[i],".csv"))
  parameter_count <- sum(sapply(motifs[[i]]$interactions, function(interaction) length(interaction$func)))

  result <- list(discrete_bifurcation = discrete_bifurcation, parameter_count = parameter_count)
  return(result)
}

#motif_qdbs <- lapply(seq_len(length(motifs)), read_result_from_file)

# Prepare the parametrisation function vectors to mark the results, to allow further manual inspection
f_1 <- sapply(pbn$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn$gene_function_vectors[[1]][[index_combination_vector[1]]], collapse =""))
})

f_2 <- sapply(pbn$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn$gene_function_vectors[[2]][[index_combination_vector[2]]], collapse =""))
})

marked_results <- list()

# Saving the results and graphs
for(i in seq_len(length(motif_qdbs))) {
  #svg(filename=paste0(motif_names[i],".svg"), width=4.8, height=4)

  print(bifurcation_plot(motif_qdbs[[i]], main = motif_nice_names[i]))

  #dev.off()

  #write.csv(motif_qdbs[[i]]$discrete_bifurcation, paste0(motif_names[i],".csv"), row.names = F)

  marked_results[[i]] <- data.frame(f_1 = f_1, f_2 = f_2,
                                    rewiring_distance = motif_qdbs[[i]]$discrete_bifurcation$rewiring_distance,
                                    attractor_landscape_similarity = motif_qdbs[[i]]$discrete_bifurcation$attractor_landscape_similarity)

  #write.csv(marked_results[[i]], paste0(motif_names[i],"_marked.csv"), row.names = F)
}

# Compute the correlations and robustness
correlation <- sapply(motif_qdbs, function(qdbf) cor(qdbf$discrete_bifurcation$rewiring_distance,
                                                     qdbf$discrete_bifurcation$attractor_landscape_similarity, method = "spearman"))

robustness <- sapply(motif_qdbs, function(qdbf) compute_rewiring_robustness(qdbf, p_elementary = 0.1))

#write.csv(data.frame(motif_names, robustness, correlation), "motifs-specific-summary.csv")
