# O(n * 2^n)
#' Create the list of all the possible binary vectors of given length, sorted lexicographically.
#' @param n A length of the vectors.
#' @return The list of the vectors.
all_binary_vectors <- function(n) {
  if (n == 0)
    return(list(numeric(0)))
  if (n == 1)
    return(list(c(0), c(1)))

  df <- expand.grid(replicate(n, 0:1, simplify = FALSE))[,n:1]
  m <- apply(df, 2, (function(x) x))
  return(unname(tapply(m,rep(1:nrow(m),ncol(m)),function(i)i)))
}

# O(n)
#' Get the number represented by binary vector (from most to least significant bit).
binary_vector_to_number <- function (binary_vector) {
  if (length(binary_vector) == 0) {
    return(0)
  }
  result <- 0
  for (bit in binary_vector) {
    result <- 2*result + bit
  }
  return(result)
}

# O(n)
#' Get the 1-based index of a binary vector representing a row in a truth-table.
binary_vector_to_one_based_index <- function (binary_vector) {
  return(binary_vector_to_number(binary_vector)+1)
}

# O(n * 2^n), n = number of nodes
#' Create the pre-to-post vector for decoding attractor states returned by Boolnet.
pre_to_post_vector <- function(states) {
  return(sapply(states, function(bv) binary_vector_to_number(rev(bv))))
}
