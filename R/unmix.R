#' @title Unmix sediment mixtures
#'
#' @description This function assesses the relative contribution of potential sediment sources to each sediment mixture in a dataset using a mass balance approach. It supports both unconstrained and constrained optimization, allowing for different methods of handling source variability.
#'
#' @param data Data frame containing sediment source and mixture data. Users should ensure their data is in a valid format by using the check_database() function before running the unmixing process.
#' @param iter The number of iterations for the variability analysis. Increase `iter` to improve the reliability and accuracy of the results. A sufficient number of iterations is reached when the output no longer changes significantly with further increases.
#' @param variability A character string specifying the type of variability to calculate. Possible values are "SD" for Standard Deviation or "SEM" for Standard Error of the Mean.
#' @param lvp A logical value to switch between classical variability analysis (lvp = FALSE) and Linear Variability Propagation (lvp = TRUE). LVP is a more accurate method for calculating uncertainty in unmixing models under high variability and extreme source apportionments.
#' @param constrained A logical value indicating whether the optimization should be constrained to physical solutions. If constrained = TRUE, the optimization will be restricted to solutions where all source contributions are within the range of 0 to 1. If constrained = FALSE, the optimization is unconstrained.
#' @param resolution An integer specifying the number of samples used in each hypercube dimension for constrained optimization. This parameter is only used when constrained = TRUE and is required to perform the analysis.
#' @param seed An integer value used to initialize the random number generator. Setting a seed ensures that the sequence of random numbers generated during the unmixing is reproducible. This is useful for debugging, testing, and comparing results across different runs. If no seed is provided, a random seed will be generated.
#'
#' @return A data frame containing the relative contributions of the sediment sources to each sediment mixture, across all iterations. The second and third rows of the result correspond to the solution for the central or mean value of the sources. The output includes an ID column to identify each mixture, a GOF (Goodness of Fit) column, and columns for each source showing their calculated contributions.
#'
#' @references
#' Latorre, B., Lizaga, I., Gaspar, L., & Navas, A. (2025). Evaluating the Impact of High Source Variability and Extreme Contributing Sources on Sediment Fingerprinting Models. *Water Resources Management*, *1-15*. https://doi.org/10.1007/s11269-025-04169-8
#' @export
unmix <- function(data, iter = 1000L, variability = "SEM", lvp = TRUE, constrained = FALSE, resolution = NA, seed = 123456L)
{
  # Check if multiple mixture samples are present in the data.
  # If so, this block processes each mixture sample individually.
	mixture_n = nrow(inputMixture(data))
	if (mixture_n > 1) {
		warning("Multiple mixture samples detected. Each sample will be processed individually.")
		 
		# Loop through each mixture sample
		for (i in 1:mixture_n) {
			# Create a subset of the data containing only one mixture
			exclude <- c()
			for (j in 1:mixture_n) {
				if(j != i) {
					exclude <- c(exclude, nrow(data)-mixture_n+j)
				}
			}
			# Subset the data to include only the current mixture sample.
			data_one_mixture <- data[-c(exclude),]
			# For the first mixture, run the unmixing model and store the initial results.
			if(i == 1) {
				results <- unmix(data_one_mixture, iter, variability, lvp, constrained, resolution, seed)
			}
			# For subsequent mixtures, run the unmixing model and append the new results to the existing ones.
			else {
				results <- rbind(results, unmix(data_one_mixture, iter, variability, lvp, constrained, resolution, seed))
			}
		}
		return(results)
	}

	if (constrained == F) {
		if(!is.na(resolution)) {
			stop("Error: Parameter 'resolution' is set but will not be used in unconstrained optimization (constrained = FALSE).")
		}
	
		if (lvp == F) {
			results <- unmix_unconstrained(data, variability = variability, iter = iter, means = is_averaged(data), seed = seed)
		} else {
			results <- unmix_unconstrained_lvp(data, variability = variability, iter = iter, means = is_averaged(data), seed = seed)
		}
		
		# Rename columns using a vector for clarity and efficiency
		colnames(results)[1:2] <- c("ID", "GOF")

		# Define the groups and mixture
		groups <- factor(data[, 2], levels = unique(data[, 2]))
		mixture <- levels(groups)[nlevels(groups)]

		# Use match() for a vectorized lookup, avoiding the slow for loop
		# This finds the row index in `data` for each ID in `results`
		row_indices <- match(mixture[results$ID], data$samples)
		results$ID <- data$ID[row_indices]
		
		# Format mixture id
		results[,1] <- paste0(mixture, ' (', results[,1], ')')
		
		# Move GOF columns to last position
		results <- results[, c(setdiff(colnames(results), "GOF"), "GOF")]
		
	  return(results)
	}
	
	if(is.na(resolution)) {
		stop("Error: Parameter 'resolution' is required for constrained optimization. A typical starting value could be resolution=20, which can then be iteratively doubled until convergence of the results.")
	}
	
  ###########################################
  sources <- inputSource(data)
  mixtures <- inputMixture(data)
  ###########################################
  
  # verify the number of sources and properties
  if (nrow(sources) -1 >= ncol(mixtures) ) {
    warning("As a minimum, n - 1 properties are required to discriminate rigorously between n sources. Additional properties are frequently required to increase the reliability of the results")
  }
  
  nsources1<- as.data.frame(unique(data[,2]))
  nsources<- nsources1[-nrow(nsources1),] 
  nsources<- as.vector(nsources) 
  
  cat("Summary of the model inputs:
    ", ncol(mixtures)-1, "variables from",nrow(nsources1)-1,"sources (",nsources,")",
      "\n")
  
  ivariability = 0
  if(variability == "SEM") {
	  ivariability = 1
  }
  
	if (lvp == F) {
		results <- unmix_c(sources, mixtures, ivariability, iter, resolution, seed)
	} else {
		results <- unmix_c_lvp(sources, mixtures, ivariability, iter, resolution, seed)
	}
  
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  # read groups (second column)
  groups <- data[, 2]
  # asume last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]
  # read sources
  sources <- data[!groups == mixture, ]
  # remove mixture level
  groups <- levels(droplevels(sources[, 2]))
  # replace column names
  colnames(results) <- c("id", "GOF", groups)
  
  # read groups (second column)
  groups <- data[, 2]
  # asume last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]
  # read mixtures
  mixtures <- data[groups == mixture, ]
  # replace sample names
  results$id <- as.character(results$id)
  # count sources rows to modify mixture id in the results
  nrow_sources <- nrow(sources)
  
  for (i in 1:ncol(results)) {
    results[, i] <- as.numeric(as.character(results[, i]))
  }
  
  results$id <- results$id + nrow_sources
  
  results <- results[order(results[, 1]), ]
  rownames(results) <- 1:nrow(results)

	cn <- colnames(results)
	cn[1] <- "ID"
	colnames(results) <- cn
	results$ID <- data[results$ID,]$ID

	# Move GOF columns to last position
	results <- results[, c(setdiff(colnames(results), "GOF"), "GOF")]
	results2 <- results[, c(setdiff(colnames(results), "GOF"))]

  if (iter==1) {  
    cat("Summary of the model outputs:",
        "\n",
        "See below the result/s of the unmixing process for the central or mean value of the sources.",
        "\n",
        "\n")
    print(aggregate(. ~ ID, data = results2, function(x) c(mean = mean(x))))
  }  
  
  else {  
    cat("Summary of the model outputs:",
        "\n",
        "The results below show the outcome of the unmixing process, which incorporated a source variability analysis \nover", iter, "iterations. Notice that the second and third rows of the result correspond to the solution for the \ncentral or mean value of the sources.",
        "\n",
        "\n")
    print(aggregate(. ~ ID, data = results2, function(x) c(mean = mean(x), SD = sd(x))))
    
  }  

	# Format mixture id
	results[,1] <- paste0(mixture, ' (', results[,1], ')')
	
  return(results)
}

