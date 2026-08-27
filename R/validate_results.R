#' @title Evaluate the mathematical consistency of a tracer selection
#'
#' @description This function assesses the mathematical consistency of a tracer selection for an apportionment result by computing the normalized error between the predicted and observed tracer concentrations in the virtual mixture. A low normalized error for all tracers indicates a consistent tracer selection. This function can be used to diagnose problems in the results of fingerprinting models.
#'
#' @param selected_data A data frame containing the characteristics of sediment sources and mixtures for the specific tracer selection to be evaluated.
#' @param apportionments A numeric vector containing the apportionment values (contributions) to be evaluated for each source, in the same order as they appear in the data.
#' @param error_threshold A numeric value (e.g., 0.05) representing the maximum acceptable 
#' normalized error. This value is used as a benchmark to categorize tracers as consistent 
#' or inconsistent in the diagnostic messages.
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
validate_results <- function(selected_data, apportionments, error_threshold = 0.05)
{
	# Input Validation
	if (!is.vector(apportionments) || is.list(apportionments)) {
		stop("Argument 'apportionments' must be a numeric vector.")
	}

	# Extract source and mixture data
	source <- inputSource(selected_data)
	mixture <- inputMixture(selected_data)
	source_n <- nrow(source)

	if (length(apportionments) != source_n) {
		stop(paste0("Apportionments length (", length(apportionments), 
		") does not match the number of sources (", source_n, ")."))
	}

	# Clean column names for internal processing
	colnames(source) <- gsub("^mean_", "", colnames(source))
	colnames(mixture) <- gsub("^mean_", "", colnames(mixture))
	
	# Error Calculation
	if(source_n == 2) {
		error_df <- CTS_error_2s(source, mixture, apportionments[1:2])
	}
	else if(source_n == 3) {
		error_df <- CTS_error_3s(source, mixture, apportionments[1:3])
	}
	else if(source_n == 4) {
		error_df <- CTS_error_4s(source, mixture, apportionments[1:4])
	}
	else if(source_n == 5) {
		error_df <- CTS_error_5s(source, mixture, apportionments[1:5])
	}
	else {
		stop(paste0("Error: Validation is not implemented for ", source_n," sources."))
	}
	
	# Scenario Analysis
	low_errors <- sum(error_df$CTS_error < error_threshold)
	total_tracers <- nrow(error_df)

	# Message logic based on the distribution of errors
	if (low_errors == total_tracers) {
		message(sprintf("The tracer selection is mathematically consistent (all normalized errors < %s); a single apportionment solution is likely present in the dataset.", 
			error_threshold))
	} 
	else if (low_errors > 0 && low_errors < total_tracers) {
		message(sprintf("Mathematical inconsistency: The apportionment solution is supported by %i tracers (normalized error < %s), but unsupported by %i tracers (normalized error >= %s). This may indicate another solution(s) present in the dataset.", 
			low_errors, error_threshold, total_tracers - low_errors, error_threshold))
	} 
	else {
		message(sprintf("Clear mathematical inconsistency: No tracers support this apportionment solution (all normalized errors >= %s). The solution may be an artifact located between multiple distinct solutions.", 
			error_threshold))
	}
	
	colnames(error_df)[colnames(error_df) == "CTS_error"] <- "normalized_error"

	return(error_df)
}

