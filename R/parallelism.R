prepare_parallel <- function() {
  num_cores <- max(parallel::detectCores() - 1, 1)
  cluster <- parallel::makePSOCKcluster(num_cores)
  doParallel::registerDoParallel(cluster, cores = num_cores)

  return(cluster)
}

parallel_apply <- function(collection, fun, packages = character(), seq_threshold = 0, log_frequency = 0) {
  log_message <- function(log_file, message) {
    write(paste(Sys.time(),message, sep = "> "), file = log_file, append = TRUE)
  }

  log_progress <- function(log_file, i, count, frequency) {
    if (frequency == 0) {
      return()
    }
    if (i %% frequency == 0) {
      log_message(log_file, paste0(100*i/count, "% done (", i, " out of ", count, ") "))
    }
  }

  execute_sequentially <- function()
  {
    if(log_frequency > 0) {
      if(is.list(collection)) {
        result <- (lapply(1:length(collection), function(i) {
          result <- fun(collection[[i]])
          log_progress(log_file, i, length(collection), log_frequency)
          return(result)
        }))
      } else {
        result <- (lapply(1:length(collection), function(i) {
          result <- fun(collection[i])
          log_progress(log_file, i, length(collection), log_frequency)
          return(result)
        }))
      }
      log_message(log_file, paste("Successfuly finished computation over", length(collection), "items."))
    } else {
      result <- (lapply(collection, fun))
    }

    return(result)
  }

  if(log_frequency > 0) {
    log_file <- paste(tempdir(), "BNRewiringRobustness.log", sep = "/")
    message(paste("Logging progress on the fly into:", log_file))
    log_message(log_file, paste("Starting computation over", length(collection), "items."))
  }

  if (seq_threshold > 0 & length(collection) < seq_threshold) {
    results <- execute_sequentially()
  } else {
    results <- tryCatch( {
        cluster <- prepare_parallel()
        `%dopar%` <- foreach::`%dopar%`

        if(log_frequency > 0) {
          if(is.list(collection)){
            result <- (
              foreach::foreach (i=1:length(collection), .packages = packages) %dopar% {
                result <- fun(collection[[i]])
                log_progress(log_file, i, length(collection), log_frequency)
                return(result)
              })
          } else {
            result <- (
              foreach::foreach (i=1:length(collection), .packages = packages) %dopar% {
                result <- fun(collection[i])
                log_progress(log_file, i, length(collection), log_frequency)
                return(result)
              })
          }
          log_message(log_file, paste("Successfuly finished computation over", length(collection), "items."))
        } else {
          result <- (foreach::foreach (i=collection, .packages = packages) %dopar% fun(i))
        }

        return(result)
      },
      error = function(cond) {
        message("Problem occured with parallel execution, restarting computation sequentially instead.")
        message("Here's the original error message:")
        message(cond)
        message()

        if(log_frequency > 0) {
          log_message(log_file, paste("Problem occured with parallel execution, restarting computation sequentially instead."))
          log_message(log_file, paste("The original error message:"))
          log_message(log_file, paste(cond))
        }

        return(execute_sequentially())
      },
      finally = {
        parallel::stopCluster(cluster)
      }
    )
  }

  return(results)
}
