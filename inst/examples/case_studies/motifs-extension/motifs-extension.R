library("BNRewiringRobustness")

#TODO test all studies files

load_network_as_RG <- function() {
  file_name <- system.file("examples", "case_studies", "motifs-extension", "network_lower_case.txt", package = "BNRewiringRobustness")
  graph <- read.csv2(file_name, header = FALSE, sep = "\t", comment.char = "#")[,c(1,3,2)]
  names(graph) <- c("n1", "edge", "n2")

  graph$edge <- sub("activator", "->", sub("repressor", "-|", sub("unknown", "-?", graph$edge)))

  return(graph)
}

#rg <- load_network_as_RG()
#write.table(rg, "e_coli.rg", row.names = FALSE, col.names = FALSE, sep = " ", quote = F)

rg <- load_RG(system.file("examples", "case_studies", "motifs-extension", "e_coli.rg", package = "BNRewiringRobustness"))

fb_vector <- sapply(1:nrow(rg), function(i) { rg$n1[i] != rg$n2[i] & any(rg$n1 == rg$n2[i] & rg$n2 == rg$n1[i])})

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

feedbacks <- data.frame(node_1, node_2, edge_1, edge_2)

pars <- rg[rg$n1 == rg$n2 & rg$edge == "->",]

crms <- feedbacks[feedbacks$node_1 %in% pars$n1 & feedbacks$node_2 %in% pars$n1
                  & feedbacks$edge_1 == "-|" & feedbacks$edge_1 == "-|",] # none

create_fb_bn <- function(i) {
  fb_rg <- data.frame(n1 = c(feedbacks$node_1[i], feedbacks$node_2[i]),
                      edge = c(feedbacks$edge_1[i], feedbacks$edge_2[i]),
                      n2 = c(feedbacks$node_2[i], feedbacks$node_1[i]))
  fb_pbn <- create_pbn_from_rg(fb_rg)

  bn <- extract_instance(fb_pbn, 1)
}

fb_bns <- lapply(seq_len(nrow(feedbacks)), create_fb_bn)

create_fb_extension <- function(i) {
  intermediate_genes <- compute_intermediate_genes(fb_bns[[i]]$genes, rg) #TODO export?
  intermediate_rg <- create_intermediate_rg(rg, intermediate_genes) #TODO export?
  # TODO actually size by edges?
  return(list(intermediate_genes = intermediate_genes, intermediate_rg = intermediate_rg, intermediate_size = nrow(intermediate_rg)))
}

fb_intermediate <- lapply(seq_len(nrow(feedbacks)), create_fb_extension)

fb_order <- order(sapply(fb_intermediate, function(fbi) fbi$intermediate_size))

compute_qdb <- function(i) {
  actual_i <- fb_order[i]
  qdb <- compute_discrete_bifurcation(fb_bns[[actual_i]], rg)
  return(qdb)
}


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

results <- list()

#TODO save time?

for (i in 3:length(fb_order)) {
  message(paste("--------- Computing: ", i, "-----------"))
  results[[i]] <- compute_qdb(i)

  write_result_to_file(i, results[[i]])
}


boxplot(results[[1]]$discrete_bifurcation$attractor_landscape_similarity~results[[1]]$discrete_bifurcation$rewiring_distance,
        data=results[[1]]$discrete_bifurcation)



read_result_from_file(1)
