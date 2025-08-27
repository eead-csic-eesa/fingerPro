#' @title Extract all possible minimal tracer combinations to identify the most discriminant.
#'
#' @description This function generates a list of all possible minimal tracer combinations
#' and serves as a crucial initial step (a "seed") in building a consistent tracer selection
#' within a sediment fingerprinting study. This analysis systematically explores various
#' minimal tracer combinations and solves the resulting determined systems of equations
#' to assess the **variability** of each combination. The **dispersion of the solution**
#' directly reflects the **discriminant capacity** of each tracer combination: a lower dispersion
#' indicates a higher discriminant capacity. While traditional methods like Discriminant
#' Function Analysis (DFA) also identify discriminant tracer combinations, this function
#' provides solutions that are **not restricted to the physically feasible space (0 < wi < 1)**.
#' This unconstrained approach is valuable for identifying problematic tracer selections
#' that might otherwise be masked when using constrained unmixing models, as discussed
#' by Latorre et al. (2021).
#'
#' @param data Data frame containing sediment source and mixtures. Users should ensure their 
#' data is in a valid format by using the check_database() function before running this function.

#' @param iter The number of iterations for the variability analysis. Increase `iter`
#'   to improve the reliability and accuracy of the results. A sufficient number of
#'   iterations is reached when the output no longer changes significantly with further increases.
#' @param seed An integer value used to initialize the random number generator.
#'   Setting a seed ensures that the sequence of random numbers generated during the unmixing
#'   is reproducible. This is useful for debugging, testing, and comparing results across
#'   different runs. If no seed is provided, a random seed will be generated.
#'
#' @return The function returns a data frame summarizing all possible tracer combinations.
#' The data frame includes the following columns for a scenario with three sources:
#' `tracers`, `w1`, `w2`, `w3`, `percent_physical`, `sd_w1`, `sd_w2`, `sd_w3`, and `max_sd_wi`.
#' Each row represents a tracer combination, detailing its corresponding solution ($w_i$),
#' the percentage of solutions that are physically feasible (0 < w_i < 1), the standard
#' deviation of the results (sd_w_i), and the maximum dispersion among all sources (max_sd_w_i).
#' The solutions are sorted in descending order, with the solution having the lowest dispersion
#' appearing first. This highlights the most discriminant combinations.
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
CTS_seeds <- function(data, iter = 1000, seed = 123456)
{
	source <- inputSource(data)
	mixture <- inputMixture(data)
	n <- nrow(source)
	
	colnames(source) <- gsub("^mean_", "", colnames(source))
	colnames(mixture) <- gsub("^mean_", "", colnames(mixture))
	
	if(n == 2) {
		seeds <- CTS_seeds_singles(source, mixture, iter, seed)
		seeds <- seeds[,c(1,2,3,6,4,5,7)]
	}
	else if(n == 3) {
		seeds <- CTS_seeds_pairs(source, mixture, iter, seed)
		seeds <- seeds[,c(1,2,3,4,8,5,6,7,9)]
	}
	else if(n == 4) {
		seeds <- CTS_seeds_triplets(source, mixture, iter, seed)
		seeds <- seeds[,c(1,2,3,4,5,10,6,7,8,9,11)]
	}
	else if(n == 5) {
		seeds <- CTS_seeds_quartets(source, mixture, iter, seed)
		seeds <- seeds[,c(1,2,3,4,5,6,12,7,8,9,10,11,13)]
	}
	else {
		stop(paste0("Error: CTS_seeds is not implemented for ", n," sources."))
	}
	
	seeds$percent_physical <- seeds$percent_physical * 100.0
	
	return(seeds)
}

