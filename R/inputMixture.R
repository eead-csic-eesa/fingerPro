#' Input sediment mixtures
#'
#' The function select and extract the sediment mixtures of the dataset.
#'
#' @param data Data frame containing source and mixtures data
#'
#' @export
#' 
inputSample <- function(data) {
  # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  mixtures <- data[data[,2] == levels(data[,2])[nlevels(data[,2])],]
  
  rownames(mixtures) <- NULL
  samples <- as.data.frame(mixtures[,-(0:2)])
  #samples <- round(samples, 2)
  
  id_samples <- c(1:nrow(samples))
  id <- paste('M',id_samples, sep='')
  
  samples <- cbind(id, samples)
  
  return(samples)
}
