library(BNRewiringRobustness)

bn <- load_BN(system.file("examples", "motifs","crm_or.bn", package = "BNRewiringRobustness"))
rg <- load_RG(system.file("examples", "rgs","general_2.rg", package = "BNRewiringRobustness"))

qdbf <- compute_discrete_bifurcation(bn, rg)
robustness <- compute_rewiring_robustness(qdbf)
print(bifurcation_plot(qdbf))
