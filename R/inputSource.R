#' Input sediment sources
#'
#' The function select and extract the source samples of the dataset.
#'
#' @param data Data frame containing source and mixtures data
#' @param na.omit Boolean to omit or not NA values when computing the mean and SD
#'
#' @export
#' 
inputSource <- function(data, na.omit = T) {

	# # If the data is averaged, extract and return only the source data
	if(is_averaged(data))
	{
		id <- head(data[,2], n=-1)
		sources1 <- data[-nrow(data), c(3:ncol(data))]
		sources <- cbind(id, sources1)
		return(sources)
	}

  # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  sources <- data[!data[,2] == levels(data[,2])[nlevels(data[,2])],]
  
  sources2 <- sources[ order(sources[,2]),]
  
  s_groups <- droplevels(sources2[,2])

  data_mean <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), mean, na.rm = na.omit)
  data_mean[[1]] <- NULL
  colnames(data_mean) <- paste('mean_', colnames(data_mean), sep='')
  
  data_sd <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), sd, na.rm = na.omit)
  
  data_sd <- as.data.frame(data_sd[,-1])
  colnames(data_sd) <- paste('sd_', colnames(data_sd), sep='')
  
  n <- data.frame(table(sources2[,2]))
  n[[1]] <- NULL
  n <- as.data.frame(n[-nrow(n),])
  colnames(n) <- 'n'
  
  id <- head(levels(data[,2]), n=-1)
  
  sources <- cbind(id, data_mean, data_sd, n)
  
  return(sources)
}
