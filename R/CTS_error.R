#' @title Evaluate the mathematical consistency of a tracer selection for an apportionment solution.
#'
#' @description This function assesses the mathematical consistency of a tracer selection for an apportionment result by computing the normalized error between the predicted and observed tracer concentrations in the virtual mixture. A low normalized error for all tracers indicates a consistent tracer selection. This function can be used to A) extend a minimal tracer combination obtained from the `CTS_seeds` function ensuring its mathematical consistency in order to select optimum tracers to perform the unmix, and to B) diagnose problems in the results of fingerprinting models.
#'
#' @param data A data frame containing the characteristics of sediment sources and mixtures.
#' @param solution A data frame or vector containing the apportionment values for each source.
#' If a data frame, the user must select a row from the CTS_seeds output based on a criteria: apportionment values should be positive (or if negative close to zero), high percentage of physically feasible solutions (percent_physical), and low dispersion indicating higher discriminant capacity.
#' If a vector, it must contain the weights for each source, in the same order as they appear in the data.
#'
#' @return A data frame containing the normalized error for each tracer.
#'
#' @details The function calculates a normalized error for each tracer to assess the consistency of a given apportionment solution. The method involves first computing a "virtual mixture" by using the proposed apportionment values to perform a weighted average of the source tracer concentrations. The error for each tracer is then the difference between the tracer concentration in the real mixture and the virtual mixture. This error is normalized by the range of the tracer, which is estimated from the extremes of the sources' confidence intervals.
#'
#' A low normalized error for all tracers (i.e., less than a predefined threshold like $0.05$) indicates a mathematically consistent tracer selection. If most tracers show low errors while a few have high errors, it suggests that those tracers may be non-conservative or less influential on the model's result. Conversely, high normalized errors in most tracers indicate mathematical inconsistency and can point to the existence of multiple partial solutions in the dataset.
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2021). A novel method for analysing
#' consistency and unravelling multiple solutions in sediment fingerprinting.
#' *Science of The Total Environment*, *789*, 147804.
#'
#' @export
CTS_error <- function(data, solution)
{
	source <- inputSource(data)
	mixture <- inputMixture(data)
	source_n <- nrow(source)
	
	colnames(source) <- gsub("^mean_", "", colnames(source))
	colnames(mixture) <- gsub("^mean_", "", colnames(mixture))
	
	if(source_n == 2 && is.data.frame(solution)) {
		return( CTS_error_2s(source, mixture, c(solution$w1, solution$w2)) )
	}
	else if(source_n == 3 && is.data.frame(solution)) {
		return( CTS_error_3s(source, mixture, c(solution$w1, solution$w2, solution$w3)) )
	}
	else if(source_n == 4 && is.data.frame(solution)) {
		return( CTS_error_4s(source, mixture, c(solution$w1, solution$w2, solution$w3, solution$w4)) )
	}
	else if(source_n == 5 && is.data.frame(solution)) {
		return( CTS_error_5s(source, mixture, c(solution$w1, solution$w2, solution$w3, solution$w4, solution$w5)) )
	}
	else if(source_n == 2 && is.vector(solution)) {
		return( CTS_error_2s(source, mixture, c(solution[1], solution[2])) )
	}
	else if(source_n == 3 && is.vector(solution)) {
		return( CTS_error_3s(source, mixture, c(solution[1], solution[2], solution[3])) )
	}
	else if(source_n == 4 && is.vector(solution)) {
		return( CTS_error_4s(source, mixture, c(solution[1], solution[2], solution[3], solution[4])) )
	}
	else if(source_n == 5 && is.vector(solution)) {
		return( CTS_error_5s(source, mixture, c(solution[1], solution[2], solution[3], solution[4], solution[5])) )
	}
	else {
		stop(paste0("Error: CTS_error is not implemented for ", source_n," sources."))
	}
}

