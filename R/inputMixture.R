#' Input sediment mixtures
#'
#' The function select and extract the sediment mixtures of the raw dataset.
#'
#' @param data Data frame containing source and mixtures data
#'
#' @export
#' 
inputMixture <- function(data) {

	# # If the data is averaged, extract and return only the mixture data
	if(is_averaged(data))
	{
		ntracer <- (ncol(data)-3)/2
		mixtures <- inputMixture(data[, c(1:(ntracer+2))])
		return(mixtures)
	}
	
  # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  mixtures <- data[data[,2] == levels(data[,2])[nlevels(data[,2])],]
  
  rownames(mixtures) <- NULL
  samples <- as.data.frame(mixtures[,-(0:2)])
  colnames(samples) <- gsub("^mean_", "", colnames(samples))
  colnames(samples) <- paste('mean_', colnames(samples), sep='')
  id <- as.vector(tail(data[,2], n=nrow(samples)))
  samples <- cbind(id, samples)
  
  return(samples)
}
