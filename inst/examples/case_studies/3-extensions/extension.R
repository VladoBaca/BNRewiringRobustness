library("BNRewiringRobustness")

######## Converting the network file to RG #########

load_network_as_RG <- function() {
  file_name <- system.file("examples", "case_studies", "3-extensions", "network_lower_case.txt", package = "BNRewiringRobustness")
  graph <- read.csv2(file_name, header = FALSE, sep = "\t", comment.char = "#")[,c(1,3,2)]
  names(graph) <- c("n1", "edge", "n2")

  graph$edge <- sub("activator", "->", sub("repressor", "-|", sub("unknown", "-?", graph$edge)))

  return(graph)
}

#rg <- load_network_as_RG()
#write.table(rg, "e_coli.rg", row.names = FALSE, col.names = FALSE, sep = " ", quote = F)

rg <- load_RG(system.file("examples", "case_studies", "3-extensions", "e_coli.rg", package = "BNRewiringRobustness"))


######## Finding the feedback loops #########

fb_vector <- sapply(1:nrow(rg), function(i) { rg$n1[i] != rg$n2[i] & any(rg$n1 == rg$n2[i] & rg$n2 == rg$n1[i])})

# the rg containing the candidate edges for FBs
rg_narrowed <- rg[fb_vector,]

node_1 <- character()
node_2 <- character()
edge_1 <- character()
edge_2 <- character()

for (i in 1:nrow(rg_narrowed)) {
  if (i < nrow(rg_narrowed) ) {
    for (j in (i+1):nrow(rg_narrowed)) {
      if(rg_narrowed$n1[i] == rg_narrowed$n2[j]
         & rg_narrowed$n1[j] == rg_narrowed$n2[i]
         & rg_narrowed$n1[i] != rg_narrowed$n2[i]
         & rg_narrowed$edge[i] != "-?"
         & rg_narrowed$edge[j] != "-?") {
        node_1 <- c(node_1, rg_narrowed$n1[i])
        node_2 <- c(node_2, rg_narrowed$n2[i])
        edge_1 <- c(edge_1, rg_narrowed$edge[i])
        edge_2 <- c(edge_2, rg_narrowed$edge[j])
      }
    }
  }
}

# the final data-frame with all the (combinations of) feedbacks
feedbacks <- data.frame(node_1, node_2, edge_1, edge_2)

# search for CRMs among the FBs
pars <- rg[rg$n1 == rg$n2 & rg$edge == "->",]

crms <- feedbacks[feedbacks$node_1 %in% pars$n1 & feedbacks$node_2 %in% pars$n1
                  & feedbacks$edge_1 == "-|" & feedbacks$edge_1 == "-|",]
# no CRMs

# coversion of fb from data-frame to a BN
create_fb_bn <- function(i) {
  fb_rg <- data.frame(n1 = c(feedbacks$node_1[i], feedbacks$node_2[i]),
                      edge = c(feedbacks$edge_1[i], feedbacks$edge_2[i]),
                      n2 = c(feedbacks$node_2[i], feedbacks$node_1[i]))
  fb_pbn <- create_pbn_from_rg(fb_rg)

  bn <- extract_instance(fb_pbn, 1)
}

fb_bns <- lapply(seq_len(nrow(feedbacks)), create_fb_bn)

######## Computation of quantitative discrete bifurcations functions #########

# computing the motif-extensions of the FB in the whole RG, to order them by computional feasibility
create_fb_extension <- function(i) {
  intermediate_genes <- compute_intermediate_genes(fb_bns[[i]]$genes, rg)
  tmp_rg <- create_intermediate_rg(rg, intermediate_genes)

  #union the edges so they are unique, making converting the dual edges into unknown
  uniq_vector <- sapply(1:nrow(tmp_rg), function(j) {
    length(which(tmp_rg$n1 == tmp_rg$n1[j] & tmp_rg$n2 == tmp_rg$n2[j])) == 1 } )

  u_n1 <- character()
  u_edge <- character()
  u_n2 <- character()

  for (j in which(!uniq_vector)) {
    if(!(any(u_n1 == tmp_rg$n1[j] & u_n2 == tmp_rg$n2[j]))) {
      u_n1 <- c(u_n1, tmp_rg$n1[j])
      u_edge <- c(u_edge, "-?")
      u_n2 <- c(u_n2, tmp_rg$n2[j])
    }
  }

  uniq <- tmp_rg[uniq_vector,]
  unioned <- data.frame(n1 = u_n1, edge = u_edge, n2 = u_n2)

  intermediate_rg <- rbind(uniq, unioned)

  return(list(intermediate_genes = intermediate_genes, intermediate_rg = intermediate_rg, intermediate_size = nrow(intermediate_rg)))
}

fb_intermediate <- lapply(seq_len(nrow(feedbacks)), create_fb_extension)

# order of FB motif extensions by number of edges in extension, from smallest
fb_order <- order(sapply(fb_intermediate, function(fbi) fbi$intermediate_size))

# computation of QDBF of the i-th smallest motif in its extension
compute_qdb <- function(i) {
  actual_i <- fb_order[i]

  qdb <- compute_discrete_bifurcation(fb_bns[[actual_i]], fb_intermediate[[actual_i]]$intermediate_rg)
  return(qdb)
}

# writing/reading the QDBF into csv and parameter count to txt (it is necessary for robustness computation)
write_result_to_file <- function(i, result) {
  write.csv(result$discrete_bifurcation, paste0("qdbf_res_",i,".csv"), quote = F, row.names = F)
  write(result$parameter_count, file=paste0("parameter_count_",i,".txt"))
}

read_result_from_file <- function(i) {
  discrete_bifurcation <- read.csv(paste0("qdbf_res_",i,".csv"))
  parameter_count <- as.numeric(readLines(paste0("parameter_count_",i,".txt")))

  result <- list(discrete_bifurcation = discrete_bifurcation, parameter_count = parameter_count)
  return(result)
}

aggregate_box <- function(discrete_bifurcation) {
  `%>%` <- magrittr::`%>%`

  return(discrete_bifurcation %>%
    ggplot2::ggplot(ggplot2::aes(y=attractor_landscape_similarity)) +
    ggplot2::geom_boxplot(fill = "gray", outlier.size = 0.05) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::ggtitle(" ") +
    ggplot2::xlab("All") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(color="white"),
                   axis.line=ggplot2::element_blank(), axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::ylab(""))
}

# the computation of QDBF, from smallest to largest motifs
results <- list()

for (i in 1:length(fb_order)) {
  message(paste("--------- Computing: ", i, "-----------"))
  results[[i]] <- compute_qdb(i)

  write_result_to_file(i, results[[i]])
}

######## Processing the results #########

# managed to compute results for smallest 5 extensions
results_count <- 5

#results <- lapply(seq_len(results_count), read_result_from_file)

# consolidating the relevant reordered data
rgs <- lapply(seq_len(results_count), function(i) fb_intermediate[[fb_order[i]]]$intermediate_rg)
gene_counts <- sapply(seq_len(results_count), function(i) length(fb_intermediate[[fb_order[i]]]$intermediate_genes))
fbs <- lapply(seq_len(results_count), function(i) feedbacks[fb_order[i],])
bns <- lapply(seq_len(results_count), function(i) fb_bns[[fb_order[i]]])

titles <- c("crp |-> fis, rg 0",
            "crp |-| fis, rg 0",
            "galr |-| gals, rg 1",
            "exur |-| uxur, rg 2",
            "rhar <-> rhas, rg 3")

# bifurcation plots
for (i in seq_len(results_count)) {
  #png(filename=paste0("extension_",i,".png"), width=4.8, height=3, units = "in", res = 300)

  print(bifurcation_plot(results[[i]], main = titles[i]))

  #dev.off()
  #png(filename=paste0("extension_aggr_",i,".png"), width=0.9, height=3, units = "in", res = 300)

  print(aggregate_box(results[[i]]$discrete_bifurcation))

  #dev.off()
}

# correlations and robustnesses
correlation <- sapply(seq_len(results_count), function(i) cor(results[[i]]$discrete_bifurcation$rewiring_distance,
                                                              results[[i]]$discrete_bifurcation$attractor_landscape_similarity,
                                                              method = "spearman"))

robustness <- sapply(seq_len(results_count), function(i) compute_rewiring_robustness(results[[i]], p_elementary = 0.1))

means <- sapply(seq_len(results_count), function(i) mean(results[[i]]$discrete_bifurcation$attractor_landscape_similarity))

#write.csv(data.frame(titles, robustness, correlation, means), "motifs-extension-summary.csv")


######## Testing different two-node networks in the first extension #########

# prepare the extension and the generation of two-node networks
rg_extension <- rgs[[1]]
rg_generator <- load_RG(system.file("examples", "rgs", "general_2.rg", package = "BNRewiringRobustness"))
pbn_generator <- create_pbn_from_rg(rg_generator)
pbn_generator$genes <- c("crp", "fis")
motifs_count <- length(pbn_generator$function_index_combinations)

motifs_robustness <- rep(0.0, motifs_count)
motifs_mean <- rep(0.0, motifs_count)


# compute the robustnesses for all two-node networks
for (i in seq_len(motifs_count)) {
  message(paste("--------- Computing: ", i, "-----------"))
  bn <- extract_instance(pbn_generator, i)

  qdb_result <- compute_discrete_bifurcation(bn, rg_extension)

  motif_robustness <- compute_rewiring_robustness(qdb_result, p_elementary = 0.1)
  motif_mean <- mean(qdb_result$discrete_bifurcation$attractor_landscape_similarity)

  motifs_robustness[i] <- motif_robustness
  motifs_mean[i] <- motif_mean
}

# the columns of the update functions vectors, to enable manual analysis of the particular BNs
f_crp <- sapply(pbn_generator$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn_generator$gene_function_vectors[[1]][[index_combination_vector[1]]], collapse =""))
})

f_fis <- sapply(pbn_generator$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn_generator$gene_function_vectors[[2]][[index_combination_vector[2]]], collapse =""))
})

#write.csv(data.frame(f_crp, f_fis, motifs_robustness, motifs_mean), "motifs-extension-comparison.csv")

#motifs_robustness <- read.csv("motifs-extension-comparison.csv")$motifs_robustness
#motifs_mean <- read.csv("motifs-extension-comparison.csv")$motifs_mean

# histogram of two-node networks robustness
draw_histogram <- function(data, line_values, line_names, main = "") {
  quantiles_values <- quantile(data, c(0.05, 0.5, 0.95))

  par(mar=c(2.3,2.3,3.8,0.5))
  h <- hist(data, breaks = 20, main = main, xlab = "", ylab = "")

  abline(v = quantiles_values, col = "blue", lty = "dashed")
  text(x = quantiles_values, y = max(h$counts)*1.1, labels = names(quantiles_values), col = "blue", mar = c(15,15,15,15), xpd = TRUE)

  draw_line <- function(x, label) {
    abline(v = x, col = "red")
    text(x = x, y = max(h$counts)*1.1, labels = label, col = "red", mar = c(15,15,15,15), xpd = TRUE, srt = 45)
  }

  for (i in seq_len(length(line_values))) {
    draw_line(line_values[i], line_names[i])
  }
}

#svg(filename = "motifs-extension-comparison.svg", width=4.8, height=3)

draw_histogram(motifs_robustness, robustness[c(1,2)], c("nfb", "dnfb"))

#dev.off()

#svg(filename = "motifs-extension-means-comparison.svg", width=4.8, height=3)

draw_histogram(motifs_mean, means[c(1,2)], c("nfb", "dnfb"))

#dev.off()


######## Testing randomized extensions #########

edges_count <- nrow(rg_extension)
genes <- create_pbn_from_rg(rg_extension)$genes

# Generate an Erdos-Renyi randomized net, so that it contains the given crp-fis feedback
generate_random_net <- function(crpfis) {
  # 1 = act, 2 = rep, 3 = unk
  regs <- c("->", "-|", "-?")
  counts <- sapply(regs, function(reg) length(which(rg_extension$edge == reg)))

  n1 <- character()
  edge <- character()
  n2 <- character()

  add_edge <- function(from, reg, to){
    n1 <<- c(n1, from)
    edge <<- c(edge, reg)
    n2 <<- c(n2, to)
    counts[reg] <<- counts[reg] - 1
  }

  #crp-fis feedback
  add_edge("fis","-|","crp")

  #the second edge of FB may be unknown, as per original GRN
  if(runif(1, 0, counts["-?"] + counts[crpfis] ) <= counts["-?"]) {
    add_edge("crp","-?","fis")
  } else {
    add_edge("crp",crpfis,"fis")
  }

  # randomly but uniquely generate all the other edges
  for(reg in regs) {
    while (counts[reg] > 0) {
      from <- sample(genes, 1)
      to <- sample(genes, 1)
      if (!any(n1 == from & n2 == to)) {
        add_edge(from, reg, to)
      }
    }
  }

  return(data.frame(n1, edge, n2))
}

# we will use a sample of 200 random extensions
samples_count <- 200

# the BNs of the two possible crp-fis FBs
bn_nfb <- bns[[1]]
bn_dnfb <- bns[[2]]

#Computation of one instance
extension_test_instance_compute <- function(i, crpfis, bn, gene_count) {
  message(paste("--------- Computing: ", i, "-----------"))

  # only consider the graphs whose extensions contain all 5 nodes
  correct_graph <- FALSE
  while(!correct_graph) {
    random_rg <- generate_random_net(crpfis)

    intermediate_genes <- compute_intermediate_genes(bn$genes, random_rg)

    correct_graph <- length(intermediate_genes) == gene_count
  }

  print(random_rg)
  qdb <- compute_discrete_bifurcation(bn, random_rg)
  rb <- compute_rewiring_robustness(qdb)
  mn <- mean(qdb$discrete_bifurcation$attractor_landscape_similarity)

  return(list(rb = rb, mn = mn))
}

# The NFB case
randomized_nfb_rb <- numeric()
randomized_nfb_mn<- numeric()

for (i in 1:samples_count) {
  res <- extension_test_instance_compute(i, "->", bn_nfb, gene_counts[1])

  randomized_nfb_rb <- c(randomized_nfb_rb, res$rb)
  randomized_nfb_mn <- c(randomized_nfb_mn, res$mn)
}

# The DNFB case
randomized_dnfb_rb <- numeric()
randomized_dnfb_mn <- numeric()

for (i in 1:samples_count) {
  res <- extension_test_instance_compute(i, "-|", bn_dnfb, gene_counts[2])

  randomized_dnfb_rb <- c(randomized_dnfb_rb, res$rb)
  randomized_dnfb_mn <- c(randomized_dnfb_mn, res$mn)
}

# Save/Load the results
#write.csv(data.frame(robustness = randomized_nfb_rb, mean = randomized_nfb_mn), "randomized-nfb-extension-comparison.csv")
#write.csv(data.frame(robustness = randomized_dnfb_rb, mean = randomized_dnfb_mn), "randomized-dnfb-extension-comparison.csv")

#randomized_nfb_rb <- read.csv("randomized-nfb-extension-comparison.csv")$robustness
#randomized_dnfb_rb <- read.csv("randomized-dnfb-extension-comparison.csv")$robustness

#randomized_nfb_mn <- read.csv("randomized-nfb-extension-comparison.csv")$mean
#randomized_dnfb_mn <- read.csv("randomized-dnfb-extension-comparison.csv")$mean

# Draw the graphs

#svg(filename = "randomized-nfb-extension-comparison.svg", width=4.8, height=3)
draw_histogram(randomized_nfb_rb, robustness[1], c("RG 0"), main = "crp |-> fis")
#dev.off()

#svg(filename = "randomized-nfb-means-extension-comparison.svg", width=4.8, height=3)
draw_histogram(randomized_nfb_mn, means[1], c("RG 0"), main = "crp |-> fis")
#dev.off()

#svg(filename = "randomized-dnfb-extension-comparison.svg", width=4.8, height=3)
draw_histogram(randomized_dnfb_rb, robustness[2], c("RG 0"), main = "crp |-| fis")
#dev.off()

#svg(filename = "randomized-dnfb-means-extension-comparison.svg", width=4.8, height=3)
draw_histogram(randomized_dnfb_mn, means[2], c("RG 0"), main = "crp |-| fis")
#dev.off()
