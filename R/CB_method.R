#' Conservative Balance (CB)
#'
#' The function read the isotopic ratio and content data of each individual tracer in your dataset and combine them together with the mixture to create scalar variables while preserving the same exact information.
#'
#' @param data Data frame containing sediment source and one mixture
#' @param Means Boolean to switch when using mean and sd data
#' @param seed Seed for the random number generator
#'
#' @return Data frame transformed to scalar variable for further analysis.
#'
#' @export
#'  
CB_Method <- function(data, Means = F, seed = 123456L)
{
  set.seed(seed)
  ###################################################################################################
  library(easypackages)
  libraries("fingerPro", "plyr", "dplyr")
  ###################################################################################################

  if (Means == F) {
    # Extract the concentrations
    df_conc0 <- data[grepl("Conc", names(data))]
    df_conc0 <- cbind(data[1:2], df_conc0)
    
    # Extract the ratios mean
    df_Rat0 <- data[!grepl("Conc", names(data))]
    
    # Remove prefixes from column names
    conc_names <- gsub("Conc_", "", names(df_conc0))
    rat_names <- gsub("R_", "", names(df_Rat0))
    
    # Check if remaining column names match
    if (identical(conc_names, rat_names)) {
      # The remaining column names match
      cat("\033[32mThe remaining column names match.\033[39m\n")
    } else {
      # The remaining column names don't match
      cat("\033[31mThe remaining column names do not match.\033[39m\n")
      
      # Find the differing column names
      differing_names <- setdiff(conc_names, rat_names)
      
      # Highlight the differing column names
      highlighted_names <- paste("\033[1m", differing_names, "\033[22m", sep = "")
      
      cat("Columns that differ: ", paste(highlighted_names, collapse = ", "), "\n")
    }
    
    # Check the number of columns
    if (ncol(df_conc0) == ncol(df_Rat0)) {
      # Same number of columns for Rat and Content
      cat("\033[32mSame number of isotopic ratio and content columns.\033[39m\n")
    } else {
      # Different number of columns for Rat and Content
      if (ncol(df_conc0) > ncol(df_Rat0)) {
        cat("\033[37mMore \033[31m\033[1misotopic ratio\033[22m\033[39m\033[37m columns than \033[31m\033[1mcontent\033[22m\033[39m\033[37m columns.\033[39m\n")
      } else {
        cat("\033[37mMore \033[31m\033[1mcontent\033[22m\033[39m\033[37m columns than \033[31m\033[1misotopic ratio\033[22m\033[39m\033[37m columns.\033[39m\n")
      }
    }
    
    # Check for NA values in df_conc0
    na_columns_conc <- names(df_conc0)[apply(is.na(df_conc0), 2, any)]
    if (length(na_columns_conc) > 0) {
      cat("\033[31mNA values detected in df_conc0 in the following column(s):\033[39m\n")
      cat(paste("\033[31m\033[1m", na_columns_conc, "\033[22m\033[39m", sep = ", "), "\n")
      cat("\033[33mPlease consider removing or replacing the NA data.\033[39m\n")
    } else {
      cat("\033[32mNo NA values detected in content columns\033[39m\n")
    }
    
    # Check for NA values in df_Rat0
    na_columns_rat <- names(df_Rat0)[apply(is.na(df_Rat0), 2, any)]
    if (length(na_columns_rat) > 0) {
      cat("\033[31mNA values detected in df_Rat0 in the following column(s):\033[39m\n")
      cat(paste("\033[31m\033[1m", na_columns_rat, "\033[22m\033[39m", sep = ", "), "\n")
      cat("\033[33mPlease consider removing or replacing the NA data.\033[39m\n")
    } else {
      cat("\033[32mNo NA values detected in isotopic ratio columns\033[39m\n")
    }
    
    df_filt <- cbind(df_Rat0, df_conc0[, c(3:ncol(df_conc0))])
    df_NT <- df_filt
    
    sr <- inputSource(df_NT)
    sm <- inputSample(df_NT)
    
    # Transform isotopes to scalar variables  
    conc <- sr[grepl("Conc", names(sr))]# Extract the concentrations
    SDconc <- conc[grepl("D", names(conc))]
    d15N <- sr[!grepl("Conc", names(sr))] # Extract the ratios mean
    d15N <- d15N[, c(2:(ncol(d15N)-1))]
    df_tmp <- cbind(d15N, conc[,c(1:(ncol(conc)/2))])
    
    # Extract the ratios from the mixture. The content of the mixture is not used as it is considered non-conservative
    mix_ratios <- sm[, c(1:((ncol(d15N)/2) + 1))] 
    
    scalar <- d15N
    for (ii in 1:(ncol(scalar)/2)) { # Check if you have even number of tracers
      
      # Calculate the mean value for each source(source ratio mean -  mixture ratios)* FA content
      scalar[ii] <- (df_tmp[,ii] - mix_ratios[,ii + 1]) * df_tmp[,ii + ncol(scalar)] 
      
      # Calculate the SD value for each source <- sqrt((SDrat * Conc_rat)^2) + ((Rat_S - Rat_Mix) * SDCon_S)^2)
      scalar[ii + ncol(scalar) / 2] <-
        sqrt(((df_tmp[, ii + ncol(scalar) / 2] * df_tmp[, ii + ncol(scalar)])^2) +
               ((df_tmp[, ii] - mix_ratios[, ii + 1]) * SDconc[,ii])^2)
    }
    
    # CB Method::Combine the transformed isotopic tracer with the other scalar variables
    scalar$Sources <- as.character(unique(df_NT[-c(nrow(df_NT)), 2])) # Create the sources column
    scalar <- scalar %>% dplyr::select("Sources", everything()) # reorder column position
    scalar$n <- sr[,ncol(sr)]
    
    # Replace numeric values with 0 in mix_ratios data frame
    scalar[nrow(scalar) + 1,] <- lapply(mix_ratios, function(x) replace(x, is.numeric(x), 0))
    
    # # Convert to NA the SD of the mixture
    scalar[nrow(scalar), c((ncol(scalar)/2 + 1):ncol(scalar))] <- NA
    
    df_unmix <- scalar
    df_unmix <- tibble::rowid_to_column(df_unmix, "ID")# Add and ID column at the beginning
    
    # Transform all to numeric just in case
    df_unmix[, 3:ncol(df_unmix)] <- as.data.frame(apply(df_unmix[, 3:ncol(df_unmix)], 2, function(x) as.numeric(as.character(x))))
    
#################################################################################################################################    
#################################################################################################################################    
  } else {
   
    # Extract the concentrations
    df_conc0 <- data[grepl("Conc", names(data))]
    df_conc0 <- cbind(data[1:2], df_conc0)
    
    # Extract the ratios mean
    df_Rat0 <- data[!grepl("Conc", names(data))]
    df_Rat0 <- df_Rat0 [, c(1:(ncol(df_Rat0)-1))]
    
    # Remove prefixes from column names
    conc_names <- gsub("Conc_", "", names(df_conc0))
    rat_names <- gsub("R_", "", names(df_Rat0))
    
    # Check if remaining column names match
    if (identical(conc_names, rat_names)) {
      # The remaining column names match
      cat("\033[32mThe remaining column names match.\033[39m\n")
    } else {
      # The remaining column names don't match
      cat("\033[31mThe remaining column names do not match.\033[39m\n")
      
      # Find the differing column names
      differing_names <- setdiff(conc_names, rat_names)
      
      # Highlight the differing column names
      highlighted_names <- paste("\033[1m", differing_names, "\033[22m", sep = "")
      
      cat("Columns that differ: ", paste(highlighted_names, collapse = ", "), "\n")
    }
    
    # Check the number of columns
    if (ncol(df_conc0) == ncol(df_Rat0)) {
      # Same number of columns for Rat and Content
      cat("\033[32mSame number of isotopic ratio and content columns.\033[39m\n")
    } else {
      # Different number of columns for Rat and Content
      if (ncol(df_conc0) > ncol(df_Rat0)) {
        cat("\033[37mMore \033[31m\033[1misotopic ratio\033[22m\033[39m\033[37m columns than \033[31m\033[1mcontent\033[22m\033[39m\033[37m columns.\033[39m\n")
      } else {
        cat("\033[37mMore \033[31m\033[1mcontent\033[22m\033[39m\033[37m columns than \033[31m\033[1misotopic ratio\033[22m\033[39m\033[37m columns.\033[39m\n")
      }
    }
    
    # Check for NA values in df_conc0
    na_columns_conc <- names(df_conc0)[apply(is.na(df_conc0[c(1:(nrow(df_conc0)-1)),]), 2, any)]
    if (length(na_columns_conc) > 0) {
      cat("\033[31mNA values detected in df_conc0 in the following column(s):\033[39m\n")
      cat(paste("\033[31m\033[1m", na_columns_conc, "\033[22m\033[39m", sep = ", "), "\n")
      cat("\033[33mPlease consider removing or replacing the NA data.\033[39m\n")
    } else {
      cat("\033[32mNo NA values detected in content columns\033[39m\n")
    }
    
    # Check for NA values in df_Rat0
    na_columns_rat <- names(df_Rat0)[apply(is.na(df_Rat0[c(1:(nrow(df_Rat0)-1)), ]), 2, any)]
    if (length(na_columns_rat) > 0) {
      cat("\033[31mNA values detected in df_Rat0 in the following column(s):\033[39m\n")
      cat(paste("\033[31m\033[1m", na_columns_rat, "\033[22m\033[39m", sep = ", "), "\n")
      cat("\033[33mPlease consider removing or replacing the NA data.\033[39m\n")
    } else {
      cat("\033[32mNo NA values detected in isotopic ratio columns\033[39m\n")
    }

    df_filt <- cbind(df_Rat0, df_conc0[, c(3:ncol(df_conc0))])
    df_NT <- df_filt
    
    sources1 <- df_NT[-nrow(df_NT), c(3:ncol(df_NT))]
    id_sources <- (1:(nrow(data) - 1))
    id <- paste('S', id_sources, sep = '')
    sr <- cbind(id, sources1)
    sm <- inputSample(df_Rat0[, c(1:((ncol(df_Rat0) / 2) + 1))])

    # Transform isotopes to scalar variables  
    conc <- sr[grepl("Conc", names(sr))]# Extract the concentrations
    SDconc <- conc[grepl("D", names(conc))]
    d15N <- sr[!grepl("Conc", names(sr))] # Extract the ratios mean
    d15N <- d15N[, c(2:ncol(d15N))]
    df_tmp <- cbind(d15N, conc[,c(1:(ncol(conc)/2))])
    
    # Extract the ratios from the mixture. The content of the mixture is not used as it is considered non-conservative
    mix_ratios <- sm 
    
    scalar <- d15N
    for (ii in 1:(ncol(scalar)/2)) { # Check if you have even number of tracers
      
      # Calculate the mean value for each source(source ratio mean -  mixture ratios)* FA content
      scalar[ii] <- (df_tmp[,ii] - mix_ratios[,ii + 1]) * df_tmp[,ii + ncol(scalar)] 
      
      # Calculate the SD value for each source <- sqrt((SDrat * Conc_rat)^2) + ((Rat_S - Rat_Mix) * SDCon_S)^2)
      scalar[ii + ncol(scalar) / 2] <-
        sqrt(((df_tmp[, ii + ncol(scalar) / 2] * df_tmp[, ii + ncol(scalar)])^2) +
               ((df_tmp[, ii] - mix_ratios[, ii + 1]) * SDconc[,ii])^2) # 42 refers to -1 from where the Con SD starts
    }
    
    scalar$Sources <- as.character(unique(df_NT[-c(nrow(df_NT)), 2])) # Create the sources column
    scalar <- scalar %>% dplyr::select("Sources", everything()) # reorder column position
    
    n_samples <- data[, ncol(data)]
    n_samples <- na.omit(as.data.frame(n_samples))
    names(n_samples) <- "n"
    scalar$n <- n_samples
    
    # Replace numeric values with 0 in mix_ratios data frame
    new_row <- c("M1", rep(0, ncol(scalar)-1))
    scalar <- rbind(scalar, new_row)
    
    # Convert to NA the SD of the mixture
    scalar[nrow(scalar), c((length(sm) + 1):ncol(scalar))] <- NA
    
    df_unmix <- scalar
    df_unmix <- tibble::rowid_to_column(df_unmix, "ID")# Add and ID column at the beginning
    
    # Transform all to numeric just in case
    df_unmix[, 3:ncol(df_unmix)] <- as.data.frame(apply(df_unmix[, 3:ncol(df_unmix)], 2, function(x) as.numeric(as.character(x))))
    
      }
  
  return(df_unmix)
  
   }
