#' Get the list of attractors and the attractor probability distributions of the states under asynchronous semantics.
get_attractors_data_async <- function(instance, instance_states, pre_to_post) {
  attractors_info <- BoolNet::getAttractors(instance, type="asynchronous", method="chosen", startStates = instance_states)

  attractors <- lapply(attractors_info$attractors, function(a) decode_attractor(a, pre_to_post))

  attractor_membership_vector <- compute_attractor_membership_vector(instance, instance_states, attractors)

  distribution_matrix <- compute_attractor_distribution_matrix(instance, instance_states, attractors, attractor_membership_vector)

  return(list(attractors = attractors, distribution_matrix = distribution_matrix))
}

#' Get the list of attractors and the attractor probability distributions of the states under synchronous semantics.
get_attractors_data_sync <- function(instance, pre_to_post) {
  attractors_info <- BoolNet::getAttractors(instance, type="synchronous", method="exhaustive", returnTable=TRUE)

  attractors <- lapply(attractors_info$attractors, function(a) decode_attractor(a, pre_to_post))

  distribution_matrix <- matrix(rep(0.0, length(pre_to_post)*length(attractors)), nrow = length(pre_to_post), ncol = length(attractors))

  for (i in 1:length(pre_to_post)) {
    sn <- pre_to_post[i]
    distribution_matrix[i,attractors_info$stateInfo$attractorAssignment[sn+1]] = 1
  }

  return(list(attractors = attractors, distribution_matrix = distribution_matrix))
}

# O(na)
#' Decode attractor from Boolnet into a vector of (binary) numbers representing states.
decode_attractor <- function(attractor, pre_to_post) {
  return(sapply(attractor$involvedStates[1,], function(as) pre_to_post[as+1]))
}

# Get the vector of states membership in attractors.
compute_attractor_membership_vector <- function(instance, instance_states, attractors) {
  attractor_membership_vector <- rep(0, length(instance_states))

  for (i in 1:length(attractors)) {
    for(state in attractors[[i]]) {
      attractor_membership_vector[state+1] = i
    }
  }

  return(attractor_membership_vector)
}

# Compute the async attractor distributions of all states in a form of matrix.
compute_attractor_distribution_matrix <- function(instance, instance_states, attractors, attractor_membership_vector) {

  # Get rid of self-debt by multiplying the other probabilities appropriately.
  self_normalize <- function(node) {
    if (debts_matrix[node, node] > 0.0) {
      node_sum <- sum(distribution_matrix[node,]) + sum(debts_matrix[node,])
      node_sum_without_self <- node_sum - debts_matrix[node, node]
      multiplier <- 1.0 / node_sum_without_self
      distribution_matrix[node,] <<- multiplier * distribution_matrix[node,]
      debts_matrix[node,] <<- multiplier * debts_matrix[node,]
      debts_matrix[node, node] <<- 0.0
    }
  }

  # The recursive node search procedure (DFS).
  process_node <- function(node) {
    # only process new nodes
    if(state_state[node] != S_NEW) {
      return()
    }

    # state is a member of attractor
    if(attractor_membership_vector[node] > 0){
      distribution_matrix[node,attractor_membership_vector[node]] <<- 1.0
      debts_matrix[node, node] <<- 0
      state_state[node] <<- S_DONE
      return()
    }

    state_state[node] <<- S_PROCESSING

    state_vector <- instance_states[[node]]

    # set self-debt to zero for now
    debts_matrix[node, node] <<- 0.0

    # iterate all successors
    for (gene in 1:n) {
      next_state <- BoolNet::stateTransition(instance, state_vector, type = "asynchronous", chosenGene = gene)

      next_state_number <- binary_vector_to_one_based_index(next_state)

      if (node != next_state_number) {
        # process recursively
        process_node(next_state_number)

        # successor is currently being processed, just add it do debts
        if(state_state[next_state_number] == S_PROCESSING){
          debts_matrix[node, next_state_number] <<- prob
        # successor is waiting on some nodes, add its prob-discounted attractor distribution and debts
        } else if(state_state[next_state_number] == S_WAITING){
          debts_matrix[node,] <<- debts_matrix[node,] + prob * debts_matrix[next_state_number,]
          distribution_matrix[node,] <<- distribution_matrix[node,] + prob * distribution_matrix[next_state_number,]
        # successor is done, add its prob-discounted attractor distribution
        } else { # S_DONE
          distribution_matrix[node,] <<- distribution_matrix[node,] + prob * distribution_matrix[next_state_number,]
        }
      } else {
        # self-loop, add to self-debt for now
        debts_matrix[node, node] <<- debts_matrix[node, node] + prob
      }
    }

    self_normalize(node)

    # if the node is waiting on some other node
    if(sum(debts_matrix[node,]) > 0.0) {
      state_state[node] <<- S_WAITING
    } else {
      state_state[node] <<- S_DONE
    }

    # for all nodes that wait on this node, pay/forward the debts
    for(node_waiting in 1:instance_states_count) {
      if (debts_matrix[node_waiting, node] > 0.0) {
        debt_prob <- debts_matrix[node_waiting, node]
        distribution_matrix[node_waiting,] <<- distribution_matrix[node_waiting,] + debt_prob * distribution_matrix[node,]
        debts_matrix[node_waiting,] <<- debts_matrix[node_waiting,] + debt_prob * debts_matrix[node,]
        debts_matrix[node_waiting, node] <<- 0.0

        if (sum(debts_matrix[node_waiting,]) == 0.0) {
          state_state[node_waiting] <<- S_DONE
        }
      }
    }

    return()
  }

  n <- length(instance$genes)
  prob <- 1.0/n
  instance_states_count <- length(instance_states)

  attractors_count <- length(attractors)

  S_NEW <- 0
  S_PROCESSING <- 1
  S_WAITING <- 2
  S_DONE <- 3

  state_state <- rep(S_NEW, instance_states_count)

  distribution_matrix <- matrix(rep(0.0,instance_states_count*attractors_count), nrow = instance_states_count, ncol = attractors_count)
  debts_matrix <- diag(instance_states_count)

  for (i in 1:instance_states_count) {
    if(state_state[i] == S_NEW){
      process_node(i)
    }
  }

  return(distribution_matrix)
}

