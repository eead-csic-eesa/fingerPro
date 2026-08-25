#' @title Compute the Conservativeness Index (CI) for individual tracers
#'
#' @description This function calculates the Conservativeness Index (CI) for each tracer based on the results of an individual tracer analysis.
#'
#' The CI index was adapted from its original definition to better describe the conservativeness of tracers in a high-dimensional space of multiple sources. The predicted source contributions from each tracer were first calculated and characterized by their centroid. Then, the CI index was calculated as the percentage of solutions with conservative apportionments (0 <= wi <= 1) relative to the centroid position. This new definition of the CI does not penalize tracers with dominant apportionments from one source and distributions close to a vertex of the physical space, unlike the previous definition.
#'
#' @param data A data frame containing the characteristics of sediment sources and mixtures.
#' @param completion_method A character string specifying the method for selecting the required remaining tracers to form a determined system of equations in the individual tracer analysis. Possible values are:
#'   "virtual": Fabricate remaining tracers virtually using generated random numbers. This method is valuable for an initial assessment of the tracer's consistency without the influence of other tracers from the dataset.
#'   "random": Randomly select remaining tracers from the dataset to complete the system. This method is useful for understanding how the tracer behaves when paired with others from the dataset.
#' @param iter The number of iterations for the variability analysis in the individual tracer analysis. Increase `iter` to improve the reliability and accuracy of the results. A sufficient number of iterations is reached when the output no longer changes significantly with further increases.
#' @param rng_init An integer value used to initialize the random number generator (RNG). Providing a starting value ensures that the sequence of random numbers generated is reproducible. This is useful for debugging, testing, and comparing results across different runs. If no value is provided, a random one will be generated.
#'
#' @return A data frame containing the CI value for each tracer.
#'
#' @references
#' Lizaga, I., Latorre, B., Bodé, S., Gaspar, L., Boeckx, P., & Navas, A. (2024). Combining isotopic and elemental tracers for enhanced sediment source partitioning in complex catchments. *Journal of Hydrology*, 631, 130768. https://doi.org/10.1016/j.jhydrol.2024.130768
#'
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2020). Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. *Science of The Total Environment*, *720*, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537
#'
#' @export
CI <- function(data, completion_method = "virtual", iter = 5000, rng_init = NULL)
{
	ita <- individual_tracer_analysis(data, completion_method, iter, rng_init)

	CI <- list()
	CI <- ita[[length(ita)]]
	names(CI)[[1]] <- "tracer"
	names(CI)[[2]] <- "CI_value"
	CI$tracer <- gsub("^mean_", "", CI$tracer)
	return(CI)
}

