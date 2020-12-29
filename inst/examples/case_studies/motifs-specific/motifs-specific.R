library("BNRewiringRobustness")

rg <- load_RG(system.file("examples", "rgs", "general_2.rg", package = "BNRewiringRobustness"))

motif_names <- c("fb_plus", "fb_minus", "fb_double_minus", "crm_and", "crm_assymetric", "crm_or")

motifs <- lapply(motif_names, function(motif_name)
  load_BN(system.file("examples", "motifs", paste0(motif_name, ".bn"), package = "BNRewiringRobustness")))

motif_qdbs <- lapply(motifs, function(motif) compute_discrete_bifurcation(motif, rg))

dev_null <- lapply(seq_len(length(motif_qdbs)), function(i) {
  #svg(filename=paste0(motif_names[i],".svg"), width=4.8, height=4)

  bifurcation_plot(motif_qdbs[[i]], poly_regression_degree = 3, main = motif_names[i])

  #dev.off()

  write.csv(motif_qdbs[[i]]$discrete_bifurcation, paste0(motif_names[i],".csv"))
})

correlation <- sapply(motif_qdbs, function(qdbf) cor(qdbf$discrete_bifurcation$rewiring_distance,
                                                      qdbf$discrete_bifurcation$attractor_landscape_similarity, method = "spearman"))

robustness <- sapply(motif_qdbs, function(qdbf) compute_rewiring_robustness(qdbf, p_elementary = 0.1))

write.csv(data.frame(motif_names, robustness, correlation), "motifs-specific-summary.csv")
