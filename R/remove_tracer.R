#' remove_tracer
#'
#' Function that excludes the properties of the sediment mixture/s selected by the user.
#'
#' @param data Data frame containing source and mixtures
#' @param rem.traz Tracers to be removed
#' 
#' @return Data frame containing sediment sources and mixtures
#'
#' @export
#'
remove_tracer <- function(data, rem.traz) {
  # remove selected tracers
  data <- data[ , !(names(data) %in% rem.traz)]
  cat(" The variable/variables that remains in your dataset is/are:", names(data[,3:ncol(data)]))
  return(data)
 }

