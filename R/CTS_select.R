#' @title Extend minimal tracer sets by evaluating mathematical consistency
#'
#' @description This function extends a minimal tracer combination obtained from the `CTS_explore` function ensuring its mathematical consistency in order to select optimum tracers to perform the unmix.
#'
#' @param data A data frame containing the characteristics of sediment sources and mixtures.
#' @param tracers_seeds A data frame containing the output from the `CTS_explore` function.
#' @param seed_id A numeric ID to select a specific row from `tracers_seeds`.
#' @param error_threshold A numeric value (e.g., 0.05). Only tracers with a normalized error below this value will be retained.
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
CTS_select <- function(data, tracers_seeds, seed_id, error_threshold = 0.05)
{
	# --- 1. Initial Checks ---
	if (is.null(tracers_seeds)) {
		stop("Argument 'tracers_seeds' is missing. Please provide the output from CTS_select.")
	}
  
	if (is.null(seed_id)) {
		stop("Argument 'seed_id' is missing. You must specify which seed to evaluate.")
	}

	# --- 2. Row Selection & Bounds Check ---
	# Check if seed_id exists in the ID column or is within row bounds
	if (!(seed_id %in% tracers_seeds$seed_id)) {
		max_id <- max(tracers_seeds$seed_id, na.rm = TRUE)
		stop(paste0("Error: seed_id '", seed_id, "' not found. ",
			"The maximum seed_id available is ", max_id, "."))
	}

	solution <- tracers_seeds[tracers_seeds$seed_id == seed_id, ]

	# Final check to ensure the resulting row isn't empty/null
	if (is.null(solution) || nrow(solution) == 0) {
		stop("The selected solution row is empty or could not be retrieved.")
	}

	source <- inputSource(data)
	mixture <- inputMixture(data)
	source_n <- nrow(source)
	
	colnames(source) <- gsub("^mean_", "", colnames(source))
	colnames(mixture) <- gsub("^mean_", "", colnames(mixture))
	
	if(source_n == 2) {
		error <- CTS_error_2s(source, mixture, c(solution$w1, solution$w2)) 
	}
	else if(source_n == 3) {
		error <- CTS_error_3s(source, mixture, c(solution$w1, solution$w2, solution$w3)) 
	}
	else if(source_n == 4) {
		error <- CTS_error_4s(source, mixture, c(solution$w1, solution$w2, solution$w3, solution$w4)) 
	}
	else if(source_n == 5) {
		error <- CTS_error_5s(source, mixture, c(solution$w1, solution$w2, solution$w3, solution$w4, solution$w5)) 
	}
	else {
		stop(paste0("Error: CTS_select is not implemented for ", source_n," sources."))
	}
	
	# --- 5. Filter Tracers based on Threshold ---
	# Identify tracers that pass the test
	consistent_tracers <- error$tracer[error$CTS_error < error_threshold]

	# Identify metadata columns
	metadata_cols <- colnames(data)[1:2]

	message(paste0("Testing seed ID ", seed_id, " (", gsub(" ", ", ", solution$tracers), ")"))
	n_selected <- length(consistent_tracers)
	if (n_selected > 0) {
		message(paste0(n_selected, " tracers have selected with a normalized error < ", error_threshold, "."))
		message(paste0("Selected tracers: ", paste(consistent_tracers, collapse = ", ")))
	} else {
		warning(paste0("No tracers met the error threshold (", error_threshold, "). Returning metadata columns only."))
		return( data[, c(metadata_cols), drop = FALSE] )
	}

	return(select_tracers(data, consistent_tracers))
}

