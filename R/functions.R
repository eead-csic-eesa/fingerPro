#' @title Check if data is averaged
#'
#' @description Checks a data frame to determine if it is formatted for averaged data. 
#' This is determined by verifying the presence of an "n" column and an equal number of columns 
#' prefixed with "mean_" and "sd_" that correspond to the same tracers.
#'
#' @param data A data frame containing sediment source and mixture data. It is expected to
#' have columns for tracer data.
#' @return A logical value. Returns `TRUE` if the data frame is formatted for averaged 
#' data (i.e., contains a column named "n" and a balanced set of "mean_" and "sd_" 
#' tracer columns). Returns `FALSE` otherwise.
#'
#' @keywords internal
is_averaged <- function(data)
{
	if (sum(grepl("^n$", tolower(colnames(data)))) == 1)
	{
		if((ncol(data)-3) %% 2 == 0)
		{
			ntracer <- (ncol(data)-3)/2
			
			names <- colnames(data)[-c(1:2)][1:(2*ntracer)]
			vector1 <- names[grepl("^mean_", names)]
			vector2 <- names[grepl("^sd_", names)]
			
			if (length(vector1) == ntracer && length(vector2) == ntracer &&
				all(gsub("^mean_", "", vector1) == gsub("^sd_", "", vector2)))
			{
				return(TRUE)
			}
		}
	}
	
	return(FALSE)
}

#' @title Check if data is formatted for isotopic analysis
#'
#' @description This function checks a data frame to determine if it is correctly formatted
#' for isotopic data, which requires a specific structure of tracer names.
#' An isotopic data frame is expected to have pairs of tracer columns,
#' where one is the raw value (or mean) and the other is a corresponding
#' content, identified by the "cont_" prefix.
#'
#' The function supports both 'raw' and 'averaged' data formats.
#'
#' @param data A data frame containing sediment source and mixture data.
#' 
#' @return A logical value. Returns `TRUE` if the data frame has the correct
#' isotopic format. Returns `FALSE` otherwise, and provides a descriptive
#' message explaining the reason for the failure.
#'
#' @keywords internal
is_isotopic <- function(data)
{
	# Get the column names from the data frame
	col_names <- colnames(data)

	# Determine if the data is averaged
	if (is_averaged(data))
	{
		#  Check for the presence of content columns
		if(!any(grepl("^mean_cont_", col_names)))
		{
			return(FALSE);
		}
		
		# Calculate the number of tracers
		num_tracers <- (length(col_names) - 3) / 2
		
		# Isotopic data requires an even number of tracers (value + content)
    if (num_tracers %% 2 == 1) {
      return(FALSE)
    }
    
    # Check if the names (without prefixes) match between the mean_ratio and mean_content columns
    if (!all(gsub("^mean_", "", col_names[3:(2 + num_tracers/2)]) == gsub("^mean_cont_", "", col_names[(3 + num_tracers/2):(2 + num_tracers)]))) {
      return(FALSE)
    }
    
    # Check if the names (without prefixes) match between the sd_ratio and sd_content columns
    if (!all(gsub("^sd_", "", col_names[(3 + num_tracers):(2 + 3*num_tracers/2)]) == gsub("^sd_cont_", "", col_names[(3 + 3*num_tracers/2):(2 + 2*num_tracers)]))) {
      return(FALSE)
    }
	}
	else
	{
		# Check for the presence of content columns
		if(!any(grepl("^cont_", col_names)))
		{
			return(FALSE);
		}
		
		# Calculate the number of tracers
    num_tracers <- (length(col_names) - 2)
    
    # Isotopic data requires an even number of tracers (value + content)
    if (num_tracers %% 2 == 1) {
      return(FALSE)
    }
    
    # Check if the tracer names (without prefixes) match between the ratio and concentration columns
    if (!all(col_names[3:(2 + num_tracers/2)] == gsub("^cont_", "", col_names[(3 + num_tracers/2):(2 + 2*num_tracers/2)]))) {
      return(FALSE)
    }
	}
	
	return(TRUE)
}

#' @title Build a raw dataset from averaged data
#'
#' @description Generates a raw (non-averaged) dataset by sampling individual observations
#' from the mean and standard deviation values provided in an averaged input data frame.
#' For each source, it generates 'n' observations for each tracer by sampling from a
#' normal distribution using the provided mean and standard deviation. Mixture data
#' is appended directly without sampling.
#'
#' @param data A data frame containing averaged source and mixture data. It is expected
#' to have columns for tracer means (prefixed with "mean_"), standard deviations
#' (prefixed with "sd_"), and a column "n" indicating the number of observations
#' for each source.
#' @return A data frame representing the raw, non-averaged dataset, with each row
#' corresponding to an individual observation.
#'
#' @export
raw_dataset <- function(data)
{
	# If the data is not averaged, just return
	if(!is_averaged(data))
	{
		return(data)
	}

	# Build a raw dataset from the averaged data
	source_n <- nrow(inputSource(data))
	tracer_n <- ncol(inputMixture(data))-1
	mixture_n <- nrow(inputMixture(data))
	
	raw_data <- data.frame(matrix(ncol = tracer_n+2, nrow = 0))
	colnames(raw_data) <- colnames(data)[1:(tracer_n+2)]
	
	k <- 1
	for (i in 1:source_n) {
		for (j in 1:inputSource(data)[i,2*(tracer_n+1)]) {
			raw_data[nrow(raw_data) + 1, ] <- NA
			raw_data[k,1] <- k
			raw_data[k,2] <- inputSource(data)[i,1]
			for (l in 1:tracer_n) {
				raw_data[k,l+2] <- rnorm(n = 1, mean = inputSource(data)[i,l+1], sd = inputSource(data)[i,l+1+tracer_n])
			}
			k <- k + 1
		}
	}
	
	for (i in 1:mixture_n) {
		raw_data[nrow(raw_data) + 1, ] <- NA
		for (j in 1:tracer_n) {
			raw_data[k,j+2] <- inputMixture(data)[i,j+1]
		}
		raw_data[k,1] <- k
		raw_data[k,2] <- inputMixture(data)[i,1]
		k <- k + 1
	}
	
	colnames(raw_data) <- gsub("^mean_", "", colnames(raw_data))
	
	return(raw_data)
}

#' @title Builds an averaged dataset from raw data
#'
#' @description Generates an averaged dataset from individual (non-averaged) observations.
#'
#' @param data A data frame containing raw source and mixture data.
#' @param na.omit Boolean to omit or not NA values when computing the mean and SD
#' @return A data frame representing the averaged dataset.
#'
#' @export
averaged_dataset <- function(data, na.omit = T)
{
	# If the data is averaged, just return
	if(is_averaged(data))
	{
		return(data)
	}

  # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  sources <- data[!data[,2] == levels(data[,2])[nlevels(data[,2])],]
  sources2 <- sources[ order(sources[,2]),]
  s_groups <- droplevels(sources2[,2])

  data_mean <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), mean, na.rm = na.omit)
  data_mean[[1]] <- NULL
  colnames(data_mean) <- paste('mean_', colnames(data_mean), sep='')
  
  data_sd <- aggregate(sources2[,3:ncol(sources2)], list(s_groups), sd, na.rm = na.omit)
  data_sd <- as.data.frame(data_sd[,-1])
  colnames(data_sd) <- paste('sd_', colnames(data_sd), sep='')
  
  n <- data.frame(table(sources2[,2]))
  n[[1]] <- NULL
  n <- as.data.frame(n[-nrow(n),])
  colnames(n) <- 'n'
  
  samples <- head(levels(data[,2]), n=-1)

  sources <- cbind(samples, data_mean, data_sd, n)
  
  mixtures <- data[data[,2] == levels(data[,2])[nlevels(data[,2])],]
  
  rownames(mixtures) <- NULL
  data_mean <- as.data.frame(mixtures[,-(0:2)])
  colnames(data_mean) <- gsub("^mean_", "", colnames(data_mean))
  colnames(data_mean) <- paste('mean_', colnames(data_mean), sep='')
  samples <- as.vector(tail(data[,2], n=nrow(data_mean)))
  data_sd <- data_mean * 0
  colnames(data_sd) <- gsub("^mean_", "sd_", colnames(data_mean))
  mixtures <- cbind(samples, data_mean, data_sd)
  mixtures$n <- 1
	
	averaged_data <- rbind(sources, mixtures)

	averaged_data$ID <- seq_len(nrow(averaged_data))
	all_cols <- colnames(averaged_data)
	averaged_data <- averaged_data[, c("ID", setdiff(all_cols, "ID"))]
	
	return(averaged_data)
}

#' @title Subset specific tracers from a dataset
#'
#' @description This function allows you to select a subset of tracer columns from a dataset.
#' It is designed to work with both isotopic and non-isotopic datasets, and also
#' with both averaged and raw data formats.
#'
#' @param data A data frame containing tracer data.
#' @param tracers A character vector of tracers to select (e.g., c("Ba", "Fe", "Cr")).
#'
#' @return A data frame containing only the specified tracer columns. The returned columns
#' will be selected based on the data format. For non-isotopic and raw data, it selects
#' the tracer columns (e.g., "tracer1"). For non-isotopic and averaged data, it selects
#' the mean and standard deviation columns (e.g., "mean_tracer1", "sd_tracer1"). For
#' isotopic and raw data, it selects the tracer and its corresponding concentration column
#' (e.g., "tracer1", "cont_tracer1"). For isotopic and averaged data, it selects the mean
#' and standard deviation for both the tracer and its concentration (e.g., "mean_tracer1",
#' "mean_cont_tracer1", "sd_tracer1", "sd_cont_tracer1").
#'
#' @keywords internal
select_tracers <- function(data, tracers) {
  names <- colnames(data)
  
  if (is_isotopic(data)) {
		if (is_averaged(data)) {
		  return(data[, c(names[1], names[2], paste0("mean_", tracers), paste0("mean_cont_", tracers), paste0("sd_", tracers), paste0("sd_cont_", tracers), names[length(names)])])
		} else {
		  return(data[, c(names[1], names[2], tracers, paste0("cont_", tracers))])
		}
  }
  else {
		if (is_averaged(data)) {
		  return(data[, c(names[1], names[2], paste0("mean_", tracers), paste0("sd_", tracers), names[length(names)])])
		} else {
		  return(data[, c(names[1], names[2], tracers)])
		}
  }
}

#' @title Read a sediment unmixing database
#'
#' @description This function automatically infers the type of sediment database
#' ("raw", "averaged", or "isotopic") based on its column names and verifies its integrity.
#' It validates column names and their order to ensure data is correctly
#' structured for subsequent package functions.
#'
#' To retain conservative tracers for subsequent analyses, it is recommended to perform
#' a minimal dataset cleaning beforehand:
#' \itemize{
#'   \item Replace BDL (below detection limit) entries with a small positive number.
#'   \item Exclude tracers whose mixture value is BDL or zero.
#'   \item Optionally, remove tracers with predominantly BDL values.
#' }
#'
#' **Database 'raw' format:**
#' This database contains individual measurements for scalar tracers. It must have the following columns
#' in order:
#' \itemize{
#'   \item{\strong{ID}:} Unique identifier for each sample.
#'   \item{\strong{samples}:} A categorical column identifying each source and mixture. The unique value representing the mixture must appear last. In cases with multiple mixture samples, they must all share the same mixture name but will be distinguished by unique entries in the \strong{ID} column.
#'   \item{\strong{tracer1, tracer2, ...}:} Columns for each tracer measurement.
#' }
#'
#' **Database 'isotopic raw' format:**
#' This database contains individual measurements for isotopic tracers, which require both
#' ratio and content data. It must have the following columns in order:
#' \itemize{
#'   \item{\strong{ID}:} Unique identifier for each sample.
#'   \item{\strong{samples}:} A categorical column identifying each source and mixture. The unique value representing the mixture must appear last. In cases with multiple mixture samples, they must all share the same mixture name but will be distinguished by unique entries in the \strong{ID} column.
#'   \item{\strong{ratio1, ratio2, ...}:} Columns with the isotopic ratio values for each tracer.
#'   \item{\strong{cont_ratio1, cont_ratio2, ...}:} Columns with the corresponding content (concentration) values for each tracer.
#' }
#'
#' **Database 'averaged' format:**
#' This database contains statistical summaries of the scalar tracer data. It must have the
#' following columns in order:
#' \itemize{
#'   \item{\strong{ID}:} Unique identifier for each sample.
#'   \item{\strong{samples}:} A categorical column identifying each source and mixture. The unique value representing the mixture must appear last. In cases with multiple mixture samples, they must all share the same mixture name but will be distinguished by unique entries in the \strong{ID} column.
#'   \item{\strong{mean_tracer1, mean_tracer2, ...}:} Columns with the mean value for each tracer.
#'   \item{\strong{sd_tracer1, sd_tracer2, ...}:} Columns with the standard deviation for each tracer.
#'   \item{\strong{n}:} The number of measurements used to calculate the mean and standard deviation.
#' }
#'
#' **Database 'isotopic averaged' format:**
#' This database contains statistical summaries for isotopic tracers. It must have the
#' following columns in order:
#' \itemize{
#'   \item{\strong{ID}:} Unique identifier for each sample.
#'   \item{\strong{samples}:} A categorical column identifying each source and mixture. The unique value representing the mixture must appear last. In cases with multiple mixture samples, they must all share the same mixture name but will be distinguished by unique entries in the \strong{ID} column.
#'   \item{\strong{mean_ratio1, mean_ratio2, ...}:} Columns with the mean isotopic ratio values.
#'   \item{\strong{mean_cont_ratio1, mean_cont_ratio2, ...}:} Columns with the mean isotopic content values.
#'   \item{\strong{sd_ratio1, sd_ratio2, ...}:} Columns with the standard deviation of the isotopic ratio values.
#'   \item{\strong{sd_cont_ratio1, sd_cont_ratio2, ...}:} Columns with the standard deviation of the isotopic content values.
#'   \item{\strong{n}:} The number of measurements.
#' }
#'
#' @param file Character string. The name of the CSV file or the path to it.
#' @param mixture Integer. The index of the mixture sample to keep if multiple are present. Defaults to 1.
#'
#' @return A data frame representing the sediment unmixing database
#'
#' @export
read_database <- function(file, mixture = 1) {

	# Determine the correct file path
	if (file.exists(file)) {
		path <- file
	} else {
		# Check if the file exists within the package's extdata folder
		path <- system.file("extdata", file, package = "fingerPro")
	}

	# Check if a valid path was actually found
	if (path == "") {
		stop("File not found locally or in package 'extdata'.")
	}

	# Read the data
	data <- read.csv(path)

	# Validate the data using your check_database function
	if (!isTRUE(check_database(data))) {
		stop()
	}
	
	# Handle multiple mixture samples
	# Get the number of rows that are identified as mixtures
	mixture_n <- nrow(inputMixture(data))

	message(paste0("Mixtures (", mixture_n, "): ", paste(tail(data[[1]], mixture_n), collapse = ", ")))
		
	if (mixture_n > 1) {
		# Check if the requested index is valid
		if (mixture > mixture_n || mixture < 1) {
			stop(paste0("Selected index (", mixture, ") is out of range. Only ", mixture_n, " mixtures found."))
		}
		
		# Calculate the starting row index of the mixture block
		total_rows <- nrow(data)
		mixture_start_idx <- total_rows - mixture_n + 1

		# Identify the rows to remove:
		# All rows from the start of the mixture block to the end...
		all_mixture_rows <- mixture_start_idx:total_rows
		# ...except for the specific one the user wants to keep
		rows_to_remove <- all_mixture_rows[-mixture]

		data <- data[-rows_to_remove, ]
		
		if(mixture == 1) {
			message(paste0("Selecting first mixture (", tail(data, 1)[,1], ") by default. Use the 'mixture' parameter to change."))
		}
		else {
			message(paste0("Selecting mixture: ", tail(data, 1)[,1]))
		}
	}

	# Find which columns (tracers) in the mixture have zero values
	# We exclude the first ID column
	mixture_data <- inputMixture(data)[, -1]
	zero_tracers <- colnames(mixture_data)[colSums(mixture_data == 0) > 0]

	if (length(zero_tracers) > 0) {
		  # Create a string of the affected tracer names
		  zero_tracers <- sub("^mean_", "", zero_tracers)
		  tracer_names <- paste(zero_tracers, collapse = ", ")
		  
		  warning(paste0("Zero values detected in the following mixture tracers: ", tracer_names, 
		                 ". Consider removing these tracers before proceeding."))
	}
	
	return(data)
}


#' @title Verify the integrity of a sediment unmixing database
#'
#' @param data A data frame to be checked.
#'
#' @return A logical value (`TRUE` if the database is valid, `FALSE` otherwise).
#'   If the check fails, the function will also print a descriptive error message.
#'
#' @keywords internal
check_database <- function(data) {
  
  # Check if the input is a data frame
  if (!is.data.frame(data)) {
    message("Error: The input 'data' is not a data frame.")
    return(FALSE)
  }

  # Check for NA values in the entire data frame
  if (any(is.na(data))) {
      na_indices <- which(is.na(data), arr.ind = TRUE)
      first_na_row <- na_indices[1, "row"]
      first_na_col <- na_indices[1, "col"]
      message(paste0("Error: The database contains missing values (NA). First NA found at row ", 
                     first_na_row, ", column ", first_na_col, " (", colnames(data)[first_na_col], ")."))
      return(FALSE)
  }

	# Check for non-numeric values at the cell level
	for (row in 1:nrow(data)) {
		for (col in 3:ncol(data)) {
		  if (!is.numeric(data[row, col])) {
		    message(paste0("Error: The database contains a non-numeric value. First non-numeric cell found at row ",
		                   row, ", column ", col, " (", colnames(data)[col], ")."))
		    return(FALSE)
		  }
		}
	}

  col_names <- colnames(data)
  
  # --- Infer database type ---
  is_averaged <- any(grepl("^mean_", col_names))
  
  # --- Check for 'averaged' database structure ---
  if (is_averaged) {
    
	  # --- Infer database type ---
	  is_isotopic <- any(grepl("^mean_cont_", col_names))
    
    # Check if the first two columns are 'ID' and 'samples'
    if (col_names[1] != "ID" || col_names[2] != "samples") {
      message("Error: In 'averaged' type, the first two columns must be 'ID' and 'samples'.")
      return(FALSE)
    }
    
    # Check if the last column is 'n'
    if (col_names[length(col_names)] != "n") {
      message("Error: In 'averaged' type, the last column must be 'n'.")
      return(FALSE)
    }
    
    # Extract tracer names (without prefixes)
    # The number of tracer-pairs is (total_columns - ID - samples - n) / 2
    num_tracers <- (length(col_names) - 3) / 2
    if (num_tracers <= 0 || num_tracers %% 1 != 0) {
      message("Error: Invalid number of columns. The 'averaged' database must have pairs of 'mean_' and 'sd_' columns.")
      return(FALSE)
    }
    
    # Check for mean_ and sd_ column names
    mean_cols <- col_names[3:(2 + num_tracers)]
    sd_cols <- col_names[(3 + num_tracers):(2 + 2 * num_tracers)]
    
    if (!all(gsub("^mean_", "", mean_cols) == gsub("^sd_", "", sd_cols))) {
      message("Error: Mismatched 'mean_' and 'sd_' tracer columns or incorrect order.")
      return(FALSE)
    }
    
    if (is_isotopic && num_tracers %% 2 == 1) {
      message("Error: Invalid number of columns. The 'averaged' isotopic database must contain mean_ratios, mean_content, sd_ratios and sd_content columns.")
      return(FALSE)
    }
    
    if (is_isotopic && !all(gsub("^mean_", "", col_names[3:(2 + num_tracers/2)]) == gsub("^mean_cont_", "", col_names[(3 + num_tracers/2):(2 + num_tracers)]))) {
      message("Error: Mismatched mean_ratios and mean_content tracer columns or incorrect order.")
      return(FALSE)
    }
    
    if (is_isotopic && !all(gsub("^sd_", "", col_names[(3 + num_tracers):(2 + 3*num_tracers/2)]) == gsub("^sd_cont_", "", col_names[(3 + 3*num_tracers/2):(2 + 2*num_tracers)]))) {
      message("Error: Mismatched sd_ratios and sd_content tracer columns or incorrect order.")
      return(FALSE)
    }
    
    # Check if the last group in 'samples' is the 'Mixture'
    if (!is.factor(data$samples)) {
        data$samples <- factor(data$samples, levels = unique(data$samples))
    }
    if (as.character(data$samples[nrow(data)]) != levels(data$samples)[nlevels(data$samples)]) {
        message("Error: The last group in the 'samples' column must represent the mixtures.")
        return(FALSE)
    }
    
    # If all checks pass
    message("The 'averaged' database is valid.")
    
    num_sources <- length(unique(data$samples))-1
    if(is_isotopic)
    {
    	message("Sources (", num_sources,"): ", paste(head(unique(data$samples), n=-1), collapse = ", "))
    	message("Tracers isotopic (", (num_tracers/2), "): ", paste(gsub("^mean_", "", col_names[3:(2 + num_tracers/2)]), collapse = ", "))
    }
    else
    {
    	message("Sources (", num_sources,"): ", paste(head(unique(data$samples), n=-1), collapse = ", "))
    	message("Tracers scalar (", num_tracers, "): ", paste(gsub("^mean_", "", mean_cols), collapse = ", "))
	  }
    
    return(TRUE)
  } 
  
  # --- Check for 'raw' database structure ---
  else {
    
	  # --- Infer database type ---
	  is_isotopic <- any(grepl("^cont_", col_names))
    
    # Check if the first two columns are 'ID' and 'samples'
    if (col_names[1] != "ID" || col_names[2] != "samples") {
      message("Error: In 'raw' type, the first two columns must be 'ID' and 'samples'.")
      return(FALSE)
    }
    
    # Check if there's at least one tracer column
    if (length(col_names) <= 2) {
      message("Error: The 'raw' database must contain at least one tracer column.")
      return(FALSE)
    }
    
    num_tracers <- (length(col_names) - 2)
    
    if (is_isotopic && num_tracers %% 2 == 1) {
      message("Error: Invalid number of columns. The 'raw' isotopic database must contain ratios and content columns.")
      return(FALSE)
    }
    
    if (is_isotopic && !all(col_names[3:(2 + num_tracers/2)] == gsub("^cont_", "", col_names[(3 + num_tracers/2):(2 + 2*num_tracers/2)]))) {
      message("Error: Mismatched ratios and content tracer columns or incorrect order.")
      return(FALSE)
    }
    
    # Check if the last group in 'samples' is the 'Mixture'
    if (!is.factor(data$samples)) {
        data$samples <- factor(data$samples, levels = unique(data$samples))
    }
    if (as.character(data$samples[nrow(data)]) != levels(data$samples)[nlevels(data$samples)]) {
        message("Error: The last group in the 'samples' column must represent the mixtures.")
        return(FALSE)
    }
    
    # If all checks pass
    message("The 'raw' database is valid.")
    
    num_sources <- length(unique(data$samples))-1
    if(is_isotopic)
    {
    	message("Sources (", num_sources,"): ", paste(head(unique(data$samples), n=-1), collapse = ", "))
    	message("Tracers isotopic (", (num_tracers/2), "): ", paste(col_names[3:(2 + num_tracers/2)], collapse = ", "))
    }
    else
    {
    	message("Sources (", num_sources,"): ", paste(head(unique(data$samples), n=-1), collapse = ", "))
    	message("Tracers scalar (", num_tracers, "): ", paste(col_names[3:(2 + num_tracers)], collapse = ", "))
	  }
    
    return(TRUE)
  }
}

