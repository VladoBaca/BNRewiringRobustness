# BNRewiringRobustness
R package for computing quantitative discrete bifurcation and rewiring robustness of Boolean networks.

For more information on usage of the package and the underlying theory, see the thesis on the subject: https://is.muni.cz/th/diu88/Thesis.pdf

Installation:

    install.packages("remotes")
    remotes::install_github("VladoBaca/BNRewiringRobustness")

Usage example:

    library(BNRewiringRobustness)
    
    bn <- load_BN(system.file("examples", "motifs","crm_or.bn", package = "BNRewiringRobustness"))
    rg <- load_RG(system.file("examples", "rgs","general_2.rg", package = "BNRewiringRobustness"))
    
    qdbf <- compute_discrete_bifurcation(bn, rg)
    robustness <- compute_rewiring_robustness(qdbf)
    print(bifurcation_plot(qdbf))
