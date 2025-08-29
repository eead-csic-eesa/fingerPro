#' @title Evaluate the mathematical consistency of a tracer selection for an apportionment solution.
#'
#' @description This function assesses the mathematical consistency of a tracer selection for an apportionment result by computing the normalized error between the predicted and observed tracer concentrations in the virtual mixture. A low normalized error for all tracers indicates a consistent tracer selection. This function can be used to diagnose problems in the results of fingerprinting models and also to extend a minimal tracer combination obtained from the `cts_error` function ensuring its mathematical consistency.
#'
#' @param source Data frame containing the sediment sources from a dataset.
#' @param mixture Data frame containing one of the dataset mixtures.
#' @param solution A vector containing the apportionment result.
#' 
#' @return Data frame containing the normalized error for each tracer.
#'
CTS_error_4s <- function(source, mixture, solution)
{
	source <- data.matrix(source[-1])
	mixture <- data.matrix(mixture[-1])

	# normalize
	cols <- (ncol(source)-1)/2
	for (col in c(1:cols))
	{
    mx <- max(source[, col] + source[, cols + col])
    mn <- min(source[, col] - source[, cols + col])
    source[, col] <- (source[, col] - mn) / (mx - mn)
    source[, cols + col] <- source[, cols + col] / (mx - mn)
    mixture[, col] <- (mixture[, col] - mn) / (mx - mn)
	}

	tracer <- colnames(source)[1:cols]
	
	CTS_error <- c()
	for (col in c(1:cols))
	{
		CTS_error[col] <- abs(source[1,col] * solution[1] + source[2,col] * solution[2] + source[3,col] * solution[3] + source[4,col] * solution[4] -  mixture[1,col])
	}
	
	return(data.frame(tracer, CTS_error))
}

