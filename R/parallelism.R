prepare_parallel <- function() {
  num_cores <- max(parallel::detectCores() - 1, 1)
  cluster <- parallel::makePSOCKcluster(num_cores)
  doParallel::registerDoParallel(cluster, cores = num_cores)

  return(cluster)
}

#TODO fix messaging?
parallel_apply <- function(collection, fun, packages = character(), message_frequency = 0, verbose = FALSE) {
  try_message <- function(i, count, frequency) {
    if (frequency == 0) {
      return()
    }
    if (i %% frequency == 0) {
      #write(paste(100*i/count, "% done (", i, " out of ", count, ") "), file = "C:\\Users\\vladimir.baca\\Documents\\test.txt", append = TRUE)
      message(paste0(100*i/count, "% done (", i, " out of ", count, ") "))
    }
  }

  results <- tryCatch( {
      cluster <- prepare_parallel()
      `%dopar%` <- foreach::`%dopar%`

      foreach_result <- foreach::foreach (i=collection, .packages = packages, .verbose = verbose) %dopar% fun(i)

      return(foreach_result)

      #if(message_frequency > 0) {
      #  if(is.list(collection)){
      #    return(foreach::foreach (i=1:length(collection), .packages = packages) %dopar% {
      #        result <- fun(collection[[i]])
      #        try_message(i, length(collection), message_frequency)
      #        return(result)
      #      })
      #  } else {
      #    return(foreach::foreach (i=1:length(collection), .packages = packages) %dopar% {
      #        result <- fun(collection[i])
      #        try_message(i, length(collection), message_frequency)
      #        return(result)
      #      })
      #  }
      #} else {
      #  return(foreach::foreach (i=collection, .packages = packages) %dopar% fun(i))
      #}
    },
    error = function(cond) {
      message("Problem occured with parallel execution, executing sequentially instead.")
      message("Here's the original error message:")
      message(cond)
      message()

      if(message_frequency > 0 & verbose) {
        if(is.list(collection)) {
          return(lapply(1:length(collection), function(i) {
            result <- fun(collection[[i]])
            try_message(i, length(collection), message_frequency)
            return(result)
          }))
        } else {
          return(lapply(1:length(collection), function(i) {
            result <- fun(collection[i])
            try_message(i, length(collection), message_frequency)
            return(result)
          }))
        }
      } else {
        return(lapply(collection, fun))
      }
    },
    finally = {
      parallel::stopCluster(cluster)
    }
  )

  return(results)
}
