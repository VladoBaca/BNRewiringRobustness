library("BNRewiringRobustness")
library("parallel")
library("doParallel")
library("foreach")
library("corrplot")

rg <- load_RG(system.file("examples", "rgs", "general_2.rg", package = "BNRewiringRobustness"))
pbn <- create_pbn_from_rg(rg)

motifs_count <- length(pbn$function_index_combinations)

compute_case <- function(case) {
  semantics <- case[[1]]
  attractor_similarity <- case[[2]]
  p_elementary <- case[[3]]

  results <- rep(0.0, motifs_count)

  for (i in seq_len(motifs_count)) {
    bn <- extract_instance(pbn, i)
    qdb_result <- compute_discrete_bifurcation(bn, rg, semantics = semantics, attractor_similarity = attractor_similarity )
    robustness <- compute_rewiring_robustness(qdb_result, p_elementary = p_elementary)

    results[i] <- robustness
  }

  return(results)
}

# the different setups of the computations
cases <- list(
  list("async", "activity", 0.001),
  list("sync", "activity", 0.001),
  list("async", "overlap", 0.001),
  list("sync", "overlap", 0.001),
  list("async", "activity", 0.1),
  list("sync", "activity", 0.1),
  list("async", "overlap", 0.1),
  list("sync", "overlap", 0.1))


num_cores <- max(parallel::detectCores() - 1, 1)
cluster <- parallel::makePSOCKcluster(num_cores)
doParallel::registerDoParallel(cluster, cores = num_cores)

results <- foreach::foreach (i=cases, .packages = c("BNRewiringRobustness")) %dopar% compute_case(i)

parallel::stopCluster(cluster)

# the columns of the update functions vectors
f_A <- sapply(pbn$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn$gene_function_vectors[[1]][[index_combination_vector[1]]], collapse =""))
})

f_B <- sapply(pbn$function_index_combinations, function(index_combination_vector) {
  return(paste(pbn$gene_function_vectors[[2]][[index_combination_vector[2]]], collapse =""))
})

# the specific motifs of interest
particular_motifs <- list(
  list("f+", c(0,1,0,1), c(0,0,1,1)),
  list("f-", c(1,0,1,0), c(0,0,1,1)),
  list("f--", c(1,0,1,0), c(1,1,0,0)),
  list("c_and", c(0,0,1,0), c(0,1,0,0)),
  list("c_asm", c(1,0,1,1), c(0,1,0,0)),
  list("c_or", c(1,0,1,1), c(1,1,0,1))
)

# the column with names of the specific motifs
motif_names <- rep("", motifs_count)

for (motif in particular_motifs) {
  row <- which(f_A == paste(motif[[2]], collapse = "") & f_B == paste(motif[[3]], collapse = ""))

  motif_names[row] <- motif[[1]]
}

# the whole data frame with all the data
results_all <- data.frame(f_A, f_B, motif_names,
                          results[[1]], results[[2]], results[[3]], results[[4]],
                          results[[5]], results[[6]], results[[7]], results[[8]])

# rename the columns with robustness values according to the setups
names(results_all)[4:11] <- sapply(cases, function(case) paste(unlist(case), collapse = "-"))

write.csv2(results_all, "motifs-2_results.csv")

# histograms of the motifs robustness for all the setups
hist_with_values_and_quantiles <- function(data, motif_names, main = "", breaks = 20, quantiles = c(0.05, 0.95)) {
  quantiles_values <- quantile(data, quantiles)
  par(mar=c(2.3,2.3,3.8,0.5))
  h <- hist(data, breaks = 20, main = main, xlab = "", ylab = "")

  abline(v = quantiles_values, col = "blue", lty = "dashed")
  text(x = quantiles_values, y = max(h$counts)*1.15, labels = names(quantiles_values), col = "blue", mar = c(15,15,15,15), xpd = TRUE)

  for (i in seq_len(length(motif_names))) {
    if (motif_names[i] != "") {

      abline(v = data[i], col = "red")
      text(x = data[i], y = max(h$counts)*1.1, labels = motif_names[i],
           col = "red", mar = c(15,15,15,15), xpd = TRUE, srt = 45)
    }
  }

  return(h)
}

for (i in seq_len(length(results))) {
  hist_with_values_and_quantiles(results[[i]], motif_names, main = paste(unlist(cases[[i]]), collapse = ", "))
}

# correlation matrix of the setups
correlation_matrix <- cor(results_all[,c(4,8,6,10,5,9,7,11)], method = "spearman")
corrplot(correlation_matrix, type = "lower", method = "number")

