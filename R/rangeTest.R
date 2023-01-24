#' Range test
#'
#' Function that excludes the properties of the sediment mixture/s outside the minimum and maximum values in the sediment sources.
#'
#' @param data Data frame containing source and mixtures
#'
#' @return Data frame containing sediment sources and mixtures
#'
#' @export
#'
rangeTest <- function(data) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  # read groups (second column)
  groups <- data[, 2]
  
  # replace all groups with 'source' except the last one
  levels(groups)[1:nlevels(groups) - 1] <- "source"
  
  # compute min values for mixture and source groups
  min <- aggregate(data[, 3:ncol(data)], list(groups), min)
  n_min <- c("s_min", "t_min")
  mint <- as.data.frame(t(min[, -1]))
  colnames(mint) <- n_min
  
  # compute max values for mixture and source groups
  max <- aggregate(data[, 3:ncol(data)], list(groups), max)
  n_max <- c("s_max", "t_max")
  maxt <- as.data.frame(t(max[, -1]))
  colnames(maxt) <- n_max
  
  # merge min and max values
  ranges1 <- merge(maxt, mint, by.x = "row.names", by.y = "row.names")
  
  # filter properties out of range
  ranges <- subset(ranges1, s_max < t_max | s_min > t_min)
  rows <- as.vector(ranges$Row.names)
  data <- data[, !(names(data) %in% rows)]
  
  cat("Attention->", nrow(ranges), "variables from a total of", nrow(ranges1), 
    "were removed:", crayon::red(crayon::bold(ranges[, 1])), ".", "\n")
  cat(" The variable/variables that remains in your dataset is/are:", 
      crayon::green(crayon::bold(names(data[, 3:ncol(data)]))), ".")
  
  return(data)
}
