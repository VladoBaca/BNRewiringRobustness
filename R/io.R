# O(n^2 * 2^(2n))
#' Load the BoolNet Boolean network from the given .bn file.
#' @param file_name File name to load from.
#' @return The Boolean network.
#' @export
load_BN <- function(file_name) {
  return(BoolNet::loadNetwork(file_name))
}

# O(2^(n * 2^n) * n)
#' Load the parametrised Boolean network from the RG in given file.
#' @param file_name File name to load from.
#' @return The parametrised Boolean network.
load_PBN <- function(file_name) {
  rg <- load_RG(file_name)
  pbn <- create_pbn_from_rg(rg)

  return(pbn)
}

#' Load the Regulatory graph from the given file.
#' @param file_name File name to load from.
#' @return The Regulatory graph.
#' @export
load_RG <- function(file_name) {
  rg_df <- read.csv2(file_name, header = FALSE, sep = " ", comment.char = "#")
  names(rg_df) <- c("n1", "edge", "n2")

  return(rg_df)
}

is_simple_string <- function(input) {
  return(is.character(input) & length(input) == 1)
}
