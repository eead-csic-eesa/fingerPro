#' @title Create a synthetic sediment mixture for validation
#'
#' @description This function generates a virtual sediment mixture based on the
#' characteristics of existing sediment sources and a set of user-defined apportionment
#' weights. It effectively simulates a mixture with known source contributions.
#'
#' @param data A data frame containing the characteristics of the sediment sources.
#' @param weights A numeric vector representing the proportional contributions (apportionment
#'   values) of each source to the virtual mixture. The order of weights in the vector must
#'   correspond to the order of sources in the `data` frame. The sum of `weights` should
#'   ideally equal 1.
#'
#' @return A data frame representing the virtual mixture. This data frame will have the
#'   same structure as a single row for a mixture in your input `data`, but with tracer
#'   values calculated based on the provided `weights`.
#'
#' @details A virtual mixture is a hypothetical sediment sample created by mathematically
#' combining the tracer characteristics of known sources according to specified proportions (`weights`).
#' This is a powerful tool in sediment fingerprinting for:
#' \itemize{
#'   \item **Consistency Checks**: Comparing observed mixture data against a virtual mixture
#'     can help assess the consistency of a dataset or the validity of an unmixing solution.
#'   \item **Scenario Testing**: Simulating mixtures under different hypothetical source contributions
#'     to understand how changes might affect sediment composition.
#'   \item **Model Validation**: Generating known virtual mixtures to test the accuracy and
#'     performance of unmixing models.
#' }
#' The function calculates the tracer values for the virtual mixture by taking the weighted
#' average of the corresponding tracer values from each source.
#'
#' @export
virtual_mixture <- function(data, weights)
{
  # Check for isotopic tracers, which are not supported by this function.
  if (any(grepl("^cont_", colnames(data)))) {
    stop("Error: This function does not support isotopic datasets.")
  }
 
  # Check if weights is a numeric vector
  if (!is.numeric(weights) || !is.vector(weights)) {
    stop("Error: 'weights' must be a numeric vector.")
  }
  
  # Check if multiple mixture samples are present in the data
  if (nrow(inputMixture(data)) > 1) {
    stop("Error: Dataset has multiple mixtures. This function works with a single mixture only.")
  }
  
	source <- inputSource(data)
	mixture <- inputMixture(data)
	source_n <- nrow(source)
	tracer_n <- (ncol(source)-2)/2

  # Check if the length of weights is equal to the number of sources
  if (length(weights) != source_n) {
    stop(paste0("Error: The length of 'weights' (", length(weights), 
		") must be equal to the number of sources (", source_n, ")."))
  }
  
  # Check if the sum of weights is approximately 1
  if (abs(sum(weights) - 1) > 1e-6) {
    warning("Warning: The sum of 'weights' does not equal 1.")
  }
  
  # Compute the tracer values for the virtual mixture
  # by calculating the weighted average of each tracer from the sources.
	for (i in 1:tracer_n) {
		avg <- 0.0
		for (j in 1:source_n) {
			avg <- avg + source[j,i+1] * weights[j]
		}
		data[nrow(data),i+2] <- avg
	}
	
  return(data)
}

