#' @title Identify minimal tracer combinations with high discriminant power
#'
#' @description This function generates a list of all possible minimal tracer combinations
#' and serves as a crucial initial step (a "seed") in building a consistent tracer selection
#' within a sediment fingerprinting study. This analysis systematically explores various
#' minimal tracer combinations and solves the resulting determined systems of equations
#' to assess the variability and reliability of each combination. The dispersion
#' of the solution directly reflects the discriminant capacity of each tracer
#' combination, where a lower dispersion indicates a higher capacity to distinguish
#' between sources. Furthermore, by evaluating solutions in an unconstrained manner,
#' the function assesses the conservativeness of the tracers; it identifies
#' whether they remain within a physically plausible range or if they exhibit
#' non-conservative behavior. While traditional methods like Discriminant Function
#' Analysis (DFA) also identify discriminant tracer combinations, this function
#' provides solutions that are not restricted to the physically feasible space (0 < wi < 1).
#' This unconstrained approach is valuable for identifying problematic tracer
#' selections that might otherwise be masked when using constrained unmixing models,
#' as discussed by Latorre et al. (2021).
#'
#' @param data Data frame containing sediment source and mixtures.
#' @param iter The number of iterations for the variability analysis. Increase `iter`
#'   to improve the reliability and accuracy of the results. A sufficient number of
#'   iterations is reached when the output no longer changes significantly with further increases.
#' @param rng_init An integer value used to initialize the random number generator (RNG). Providing a starting value ensures that the sequence of random numbers generated is reproducible. This is useful for debugging, testing, and comparing results across different runs. If no value is provided, a random one will be generated.
#'
#' @return The function returns a data frame summarizing all possible tracer combinations.
#' The data frame includes the following columns for a scenario with three sources:
#' `tracers`, `w1`, `w2`, `w3`, `percent_physical`, `sd_w1`, `sd_w2`, `sd_w3`, and `max_sd_wi`.
#' Each row represents a tracer combination, detailing its corresponding solution ($w_i$),
#' the percentage of solutions that are physically feasible (0 < w_i < 1), the standard
#' deviation of the results (sd_w_i), and the maximum dispersion among all sources (max_sd_w_i).
#' The solutions are sorted in descending order, with the solution having the lowest dispersion
#' appearing first. This highlights the most discriminant and conservative combinations.
#'
#' @details The Consistent Tracer Selection (CTS) method, as described by Latorre et al. (2021),
#' begins by considering all possible sets of $n-1$ tracers, where $n$ is the number of sources.
#' Each of these sets forms a determined system of linear equations that can be solved.
#' To account for the variability within the sources, each tracer set is iteratively solved.
#' This process involves sampling the source average values from a t-distribution, reflecting
#' the discrepancy between the true mean and the measured mean due to finite observations.
#' The maximum dispersion observed in the average apportionments for each tracer set is then
#' used as a criterion to rank them, with lower dispersion indicating higher discriminant capacity.
#' This initial step is crucial for identifying multiple discriminant solutions within the dataset,
#' a problem often unexplored by traditional tracer selection methods.
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2021). A novel method for analysing
#' consistency and unravelling multiple solutions in sediment fingerprinting.
#' *Science of The Total Environment*, *789*, 147804.
#'
#' @export
CTS_explore <- function(data, iter = 1000, rng_init = NULL)
{
	source <- inputSource(data)
	mixture <- inputMixture(data)
	source_n <- nrow(source)
	tracer_n <- ncol(mixture)-1
	
	colnames(source) <- gsub("^mean_", "", colnames(source))
	colnames(mixture) <- gsub("^mean_", "", colnames(mixture))

	# Minimal combinations use (source_n - 1) tracers
	k <- source_n - 1
	total_combos <- choose(tracer_n, k)

	# Text labels for the message
	combo_label <- switch(as.character(k), 
		"1" = "singles", 
		"2" = "pairs", 
		"3" = "triplets", 
		"4" = "quartets", 
		"combinations")

	message(sprintf("Computing all tracer %s, %i in total.", combo_label, total_combos))

	if(source_n == 2) {
		tracers_seeds <- CTS_seeds_singles(source, mixture, iter, rng_init)
		tracers_seeds <- tracers_seeds[,c(1,2,3,6,4,5,7)]
	}
	else if(source_n == 3) {
		tracers_seeds <- CTS_seeds_pairs(source, mixture, iter, rng_init)
		tracers_seeds <- tracers_seeds[,c(1,2,3,4,8,5,6,7,9)]
	}
	else if(source_n == 4) {
		tracers_seeds <- CTS_seeds_triplets(source, mixture, iter, rng_init)
		tracers_seeds <- tracers_seeds[,c(1,2,3,4,5,10,6,7,8,9,11)]
	}
	else if(source_n == 5) {
		tracers_seeds <- CTS_seeds_quartets(source, mixture, iter, rng_init)
		tracers_seeds <- tracers_seeds[,c(1,2,3,4,5,6,12,7,8,9,10,11,13)]
	}
	else {
		stop(paste0("Error: CTS_explore is not implemented for ", n," sources."))
	}
	
	tracers_seeds$percent_physical <- tracers_seeds$percent_physical * 100.0

	# Create the seed_id column
	tracers_seeds$seed_id <- seq_len(nrow(tracers_seeds))
	all_cols <- colnames(tracers_seeds)
	tracers_seeds <- tracers_seeds[, c("seed_id", setdiff(all_cols, "seed_id"))]
	
	return(tracers_seeds)
}

