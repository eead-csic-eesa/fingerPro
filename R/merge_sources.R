#' merge_sources
#'
#' Function that merge two sources into one.
#'
#' @param data Data frame containing source and mixtures
#' @param S0 First source to be merged
#' @param S1 Second source to be merged
#' @return Data frame containing sediment sources and mixtures
#'
#' @export
#'

merge_sources <- function(data, S0, S1) {
  # merge both selected sources S0 + S1
  data[, 2][data[, 2] == S0] <- S1
  return(data)
}