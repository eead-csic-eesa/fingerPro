#' Conservative Balance (CB) Method for Isotopic Tracer Analysis
#'
#' This function transforms isotopic ratio and content data of individual tracers
#' in a dataset into virtual elemental tracers, which can then be combined
#' with classical tracers and analyzed with standard unmixing models.
#'
#' @param data A data frame containing the isotopic tracer characteristics of
#' sediment sources and mixtures. The data should be correctly formatted for
#' isotopic analysis, including both isotopic ratio and isotopic content. Users
#' should ensure their data is in a valid format by using the check_database() 
#' function before running the CB method.
#'
#' @return A data frame where isotopic tracers have been converted into
#' scalar virtual tracers for further analysis. After the transformation,
#' the mixture's row will have tracer values of zero.
#'
#' @details The Conservative Balance (CB) method provides a novel,
#' physically-based framework for analyzing isotopic tracers in sediment
#' fingerprinting.
#'
#' The core of the method is an exact transformation that combines the isotopic
#' ratio and isotopic content into a virtual elemental tracer. This approach
#' has two key advantages: it allows isotopic tracers to be analyzed using
#' classical unmixing models, and it enables their combined use with elemental
#' tracers to potentially increase the discriminant capacity of the
#' fingerprinting analysis.
#'
#' This function implements the simplified approximation of the CB transformation,
#' assuming that the isotopic ratio is much smaller than 1.
#' The calculation is performed for both averaged and non-averaged datasets.
#'
#' A key feature of this transformation is that the tracer values for the
#' mixture are set to zero. This is a direct consequence of the method, as the
#' isotopic ratio of each source is subtracted from the mixture's isotopic ratio,
#' meaning the mixture's own value minus itself results in zero.
#'
#' @references
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2022). Combined use of 
#  geochemistry and compound-specific stable isotopes for sediment
#' fingerprinting and tracing. Science of The Total Environment, 832, 154834.
#'
#' @export
CB_method <- function(data)
{
	# Check format
	if (!is_isotopic(data))
	{
		stop("The dataset is not correctly formatted for isotopic analysis. Please check the column names and structure.")
	}

  # Check for NA values in the entire data frame
  if (any(is.na(data))) {
    na_indices <- which(is.na(data), arr.ind = TRUE)
    first_na_row <- na_indices[1, "row"]
    first_na_col <- na_indices[1, "col"]
    stop(paste0("Error: The database contains missing values (NA). First NA found at row ", 
     first_na_row, ", column ", first_na_col, " (", colnames(data)[first_na_col], ")."))
  }

  if(nrow(inputMixture(data))==0)
  {
  	stop("The dataset does not contain a mixture. The mixture's isotopic ratios are required for the CB method.")
  }
  
  if(nrow(inputMixture(data))>1)
  {
  	stop("The dataset contains multiple mixtures. Only one mixture is supported by the CB method.")
  }

	# Get the column names from the data frame
	col_names <- colnames(data)
	
	# Determine if the data is averaged
  if (is_averaged(data))
  {
		# Calculate the number of tracers
		num_tracers <- (length(col_names) - 3) / 4
		
    # Extract the mean source ratios
    source_mean_ratio <- data[,c(1,2,3:(2+num_tracers))]

    # Extract the mean source content
		source_mean_content <- data[,c(1,2,(3+num_tracers):(2+2*num_tracers))]

    # Extract the sd source ratios
    source_sd_ratio <- data[,c(1,2,(3+2*num_tracers):(2+3*num_tracers))]

    # Extract the sd source content
		source_sd_content <- data[,c(1,2,(3+3*num_tracers):(2+4*num_tracers))]

    # Extract the mixture ratios
		mix_ratio <- data[nrow(data),][,c(1,2,3:(2+num_tracers))]
		
		# Create the CB data frame		
		cb <- data[,c(1,2,3:(2+num_tracers),(3 + 2*num_tracers):(2 + 3*num_tracers),ncol(data))]
		colnames(cb)[3:(2+num_tracers)] <- gsub("^mean_", "mean_CB_", colnames(cb)[3:(2+num_tracers)])
		colnames(cb)[(3+num_tracers):(2+2*num_tracers)] <- gsub("^sd_", "sd_CB_", colnames(cb)[(3+num_tracers):(2+2*num_tracers)])
		
		# Calculate the equivalent scalar tracer for each source ratio and content pair, by subtracting the mixture's ratio
		for (i in 1:nrow(cb)) {
			for (j in 3:(2+num_tracers)) {
				cb[i,j] <- source_mean_content[i,j] * (source_mean_ratio[i,j] - mix_ratio[1,j])
			}
		}
		
		# Calculate the variability of the equivalent scalar tracer
		for (i in 1:nrow(cb)) {
			for (j in 3:(2+num_tracers)) {
				cb[i,j+num_tracers] <- sqrt( ( source_sd_content[i,j] * (source_mean_ratio[i,j] - mix_ratio[1,j]) )^2 + 
					( source_mean_content[i,j] * source_sd_ratio[i,j] )^2 )
			}
		}
         
    return(cb)
  }
  else
  {
		# Calculate the number of tracers
		num_tracers <- (length(col_names) - 2) / 2

    # Extract the source ratios
    source_ratio <- data[,c(1,2,3:(2+num_tracers))]

    # Extract the source content
		source_content <- data[,c(1,2,(3+num_tracers):(2+2*num_tracers))]

    # Extract the mixture ratios
		mix_ratio <- data[nrow(data),][,c(1,2,3:(2+num_tracers))]

		# Create the CB data frame		
		cb <- source_ratio
		colnames(cb) <- gsub("^", "CB_", colnames(cb))
		colnames(cb)[1] <- colnames(source_ratio)[1]
		colnames(cb)[2] <- colnames(source_ratio)[2]
		
		# Calculate the equivalent scalar tracer for each source ratio and content pair, by subtracting the mixture's ratio
		for (i in 1:nrow(cb)) {
			for (j in 3:ncol(cb)) {
				cb[i,j] <- source_content[i,j] * (source_ratio[i,j] - mix_ratio[1,j])
			}
		}
		
		return(cb)
  }
}
