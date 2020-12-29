library("BNRewiringRobustness")

rg <- load_RG(system.file("examples", "rgs", "general_2.rg", package = "BNRewiringRobustness"))

motif_names <- c("fb_plus", "fb_minus", "fb_double_minus", "crm_and", "crm_assymetric", "crm_or")

motifs <- lapply(motif_names, function(motif_name)
  load_BN(system.file("examples", "motifs", paste0(motif_name, ".bn"), package = "BNRewiringRobustness")))

motif_qdbs <- lapply(motifs, function(motif) compute_discrete_bifurcation(motif, rg))

read_result_from_file <- function(i) {
  discrete_bifurcation <- read.csv(paste0(motif_names[i],".csv"))
  parameter_count <- sum(sapply(motifs[[i]]$interactions, function(interaction) length(interaction$func)))

  result <- list(discrete_bifurcation = discrete_bifurcation, parameter_count = parameter_count)
  return(result)
}

#motif_qdbs <- lapply(seq_len(length(motifs)), read_result_from_file)

for(i in seq_len(length(motif_qdbs))) {
  #svg(filename=paste0(motif_names[i],".svg"), width=4.8, height=4)

  print(bifurcation_plot(motif_qdbs[[i]], main = motif_names[i]))

  #dev.off()

  #write.csv(motif_qdbs[[i]]$discrete_bifurcation, paste0(motif_names[i],".csv"), row.names = F)
}

correlation <- sapply(motif_qdbs, function(qdbf) cor(qdbf$discrete_bifurcation$rewiring_distance,
                                                      qdbf$discrete_bifurcation$attractor_landscape_similarity, method = "spearman"))

robustness <- sapply(motif_qdbs, function(qdbf) compute_rewiring_robustness(qdbf, p_elementary = 0.1))

write.csv(data.frame(motif_names, robustness, correlation), "motifs-specific-summary.csv")
