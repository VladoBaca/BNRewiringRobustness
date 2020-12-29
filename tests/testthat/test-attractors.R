context("Functionality of constructing attractor probability distributions")

getAsyncAttractorsData_iteratively <- function(network, attractors) {
  stateInfo <- computeAttractorProbabilities_iteratively(network, attractors)

  return(stateInfo)
}

computeAttractorProbabilities_iteratively <- function(network, attractors) {
  n <- length(network$genes)

  allStates <- all_binary_vectors(n)
  allStatesCount <- 2^n

  attractorsCount <- length(attractors)

  attractorMembershipVector <- rep(0, allStatesCount)

  for (i in 1:attractorsCount) {
    for(state in attractors[[i]]) {
      attractorMembershipVector[state+1] = i
    }
  }

  pureStates <- which(sapply(0:(allStatesCount-1), function(s) attractorMembershipVector[s+1] == 0)) - 1
  pureStatesCount <- length(pureStates)

  pureStateIdentityVector <- rep(0, allStatesCount)

  if (pureStatesCount > 0) {
    for (i in 1:pureStatesCount) {
      pureState <- pureStates[i]
      pureStateIdentityVector[pureState+1] <- i
    }
  }

  matrixDim <- pureStatesCount + attractorsCount
  mat <- matrix(rep(0.0,matrixDim*matrixDim), nrow = matrixDim, ncol = matrixDim)

  for (pureState in pureStates) {
    stateVector <- allStates[[pureState+1]]
    successors <- numeric(0)

    for (i in 1:n) {
      nextState <- BoolNet::stateTransition(network, stateVector, type = "asynchronous", chosenGene = i)
      if (!all(stateVector == nextState)) {
        nextStateNumber <- binary_vector_to_number(nextState)

        if (attractorMembershipVector[nextStateNumber+1] != 0) {
          successors <- append(successors, pureStatesCount + attractorMembershipVector[nextStateNumber+1])
        }
        else {
          successors <- append(successors, pureStateIdentityVector[nextStateNumber+1])
        }
      }
    }

    successorsCount <- length(successors)
    unitProbability <- 1.0 / successorsCount

    pureStateIndex <- pureStateIdentityVector[pureState+1]

    for (s in successors) {
      mat[pureStateIndex, s] <- mat[pureStateIndex, s] + unitProbability
    }
  }

  for (i in 1:attractorsCount) {
    mat[pureStatesCount + i,pureStatesCount + i] <- 1.0
  }

  epsilon <- 0.00000001
  minimalValue <- 1 - epsilon

  checkConvergence <- function() {
    if (pureStatesCount > 0){
      for (i in 1:pureStatesCount) {
        if (sum(mat[i,(pureStatesCount+1):matrixDim]) <= minimalValue) {
          return(FALSE)
        }
      }
    }
    return(TRUE)
  }

  i <- 0
  while (!checkConvergence()) {
    mat <- mat %*% mat
    i <- i+1
  }

  stateInfo <- lapply(1:allStatesCount, function(s) {
    if (attractorMembershipVector[s] != 0) {
      return(list(attractors = c(attractorMembershipVector[s]), probabilities = c(1.0)))
    }
    else {
      attractorProbs <- mat[pureStateIdentityVector[s],(pureStatesCount+1):matrixDim]
      boolVector <- attractorProbs != 0
      return(list(attractors = which(boolVector), probabilities = attractorProbs[boolVector]))
    }
  })

  return(stateInfo)
}

# Transforms attractor distribution matrix to list of attractors and probabilites of each state
attractor_distribution_matrix_to_sparse_list <- function(distribution_matrix) {

  sparse_list <- lapply(1:nrow(distribution_matrix), function(state) {
    attractor_probs <- distribution_matrix[state,]
    attractor_present_bool_vector <- attractor_probs != 0
    return(list(attractors = which(attractor_present_bool_vector), probabilities = attractor_probs[attractor_present_bool_vector]))
  })

  return(sparse_list)
}

test_instance <- function(instance_name) {
  instance <- load_instance_by_name(instance_name)

  dfs_attractors_data <- get_dfs_attractors_data(instance)

  dfs_result <- attractor_distribution_matrix_to_sparse_list(dfs_attractors_data$distribution_matrix)
  iterative_result <- getAsyncAttractorsData_iteratively(instance, dfs_attractors_data$attractors)

  expect_equivalent(dfs_result, iterative_result)
}

load_instance_by_name <- function(instance_name) {
  return(load_BN(system.file("examples", "motifs", paste0(instance_name, ".bn"), package = "BNRewiringRobustness")))
}

get_dfs_attractors_data <- function(instance) {
  instance_gene_count <- length(instance$genes)
  instance_states <- all_binary_vectors(instance_gene_count)
  pre_to_post <- pre_to_post_vector(instance_states)

  dfs_attractors_data <- get_attractors_data_async(instance, instance_states, pre_to_post)

  return(dfs_attractors_data)
}

get_dfs_attractors_data_sync <- function(instance) {
  instance_gene_count <- length(instance$genes)
  instance_states <- all_binary_vectors(instance_gene_count)
  pre_to_post <- pre_to_post_vector(instance_states)

  dfs_attractors_data <- get_attractors_data_sync(instance, pre_to_post)

  return(dfs_attractors_data)
}

test_that("get_attractors_data_async creates the same attractor distributions as the iterative algorithm", {
  test_instance("eq_test")
  test_instance("fb_minus")
  test_instance("ffl_c1")
  test_instance("crm_assymetric")
  test_instance("cellcycle_and")
  test_instance("mixed_regulated_double_minus")
  test_instance("fb_minus")
  test_instance("cellcycle_or")
  test_instance("interesting")
  test_instance("complex_test")
})

test_that("get_attractors_data_async computes correct attractors and distributions", {
  dfs_attractors_data <- get_dfs_attractors_data(load_instance_by_name("two_attractors"))

  unlisted_attractors <- unlist(dfs_attractors_data$attractors)
  sorted_attractors <- sort(unlisted_attractors)

  expect_equal(length(dfs_attractors_data$attractors), 2)
  expect_equivalent(sorted_attractors, c(3, 6))

  encoder <- match(sorted_attractors, unlisted_attractors)

  expected_attractors <- list(c(encoder[2]), c(1,2), c(1,2), c(encoder[1]), c(encoder[2]), c(1,2), c(encoder[2]), c(1,2))
  expected_probs <- list(c(1), c(5/12,7/12)[encoder], c(1/3,2/3)[encoder], c(1), c(1), c(1/4,3/4)[encoder], c(1), c(1/2,1/2)[encoder])

  actual_sparse_list <- attractor_distribution_matrix_to_sparse_list(dfs_attractors_data$distribution_matrix)

  actual_attractors <- lapply(actual_sparse_list, function(s) s$attractors)
  actual_probs <- lapply(actual_sparse_list, function(s) s$probabilities)

  expect_equivalent(actual_attractors, expected_attractors)
  expect_equivalent(actual_probs, expected_probs)
})

test_that("get_attractors_data_async computes correct attractors and distributions for test_triple_or", {
  dfs_attractors_data <- get_dfs_attractors_data(load_instance_by_name("test_triple_or"))

  expected_attractors <- list(c(7),c(0,2))
  expected_matrix <- matrix(c(0, 1,
                              0.5,0.5,
                              0, 1,
                              0.5,0.5,
                              0.5,0.5,
                              1, 0,
                              0.5,0.5,
                              1, 0), nrow = 8, ncol = 2, byrow = TRUE)

  if(dfs_attractors_data$attractors[[1]] != c(7)) {
    expected_attractors <- list(c(0,2),c(7))
    expected_matrix <- expected_matrix[,c(2,1)]
  }

  expect_equivalent(length(dfs_attractors_data$attractors), length(expected_attractors))
  expect_equivalent(sort(dfs_attractors_data$attractors[[1]]), expected_attractors[[1]])
  expect_equivalent(sort(dfs_attractors_data$attractors[[2]]), expected_attractors[[2]])

  expect_equivalent(dfs_attractors_data$distribution_matrix, expected_matrix)
})

test_that("get_attractors_data_sync computes correct attractors and distributions for test_triple_or", {
  dfs_attractors_data <- get_dfs_attractors_data_sync(load_instance_by_name("test_triple_or"))

  expected_attractors <- list(c(7),c(0,2), c(3,4))

  expected_matrix <- matrix(c(0, 1, 0,
                              0, 0, 1,
                              0, 1, 0,
                              0, 0, 1,
                              0, 0, 1,
                              1, 0, 0,
                              0, 0, 1,
                              1, 0, 0
                              ), nrow = 8, ncol = 3, byrow = TRUE)

  encoding_vector <- sapply(seq_len(length(expected_attractors)), function(i) {
    for (j in seq_len(length(dfs_attractors_data$attractors))) {
      if(all(sort(dfs_attractors_data$attractors[[j]]) == expected_attractors[[i]])) {
        return(j)
      }
    }
  })

  expected_attractors <- list(expected_attractors[[encoding_vector[1]]],
                              expected_attractors[[encoding_vector[2]]],
                              expected_attractors[[encoding_vector[3]]])
  expected_matrix <- expected_matrix[,encoding_vector]

  expect_equivalent(length(dfs_attractors_data$attractors), length(expected_attractors))
  expect_equivalent(sort(dfs_attractors_data$attractors[[1]]), expected_attractors[[1]])
  expect_equivalent(sort(dfs_attractors_data$attractors[[2]]), expected_attractors[[2]])
  expect_equivalent(sort(dfs_attractors_data$attractors[[3]]), expected_attractors[[3]])

  expect_equivalent(dfs_attractors_data$distribution_matrix, expected_matrix)
})
