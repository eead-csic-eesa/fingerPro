#' @title Individual tracer analysis
#'
#' @description This function computes the distribution of apportionments compatible with each individual tracer in the dataset, providing insights into the tracer's discriminant capacity and conservativeness. The method assesses the contribution of a single tracer to an unmixing model by solving a determined system of equations for each tracer.
#'
#' @param data A data frame containing the characteristics of sediment sources and mixtures. Users should ensure their data is in a valid format by using the `check_database()` function before running the individual tracer analysis.
#' @param completion_method A character string specifying the method for selecting the required remaining tracers to form a determined system of equations. Possible values are:
#'   "virtual": Fabricate remaining tracers virtually using generated random numbers. This method is valuable for an initial assessment of the tracer's consistency without the influence of other tracers from the dataset.
#'   "random": Randomly select remaining tracers from the dataset to complete the system. This method is useful for understanding how the tracer behaves when paired with others from the dataset.
#' @param iter The number of iterations for the variability analysis. Increase `iter` to improve the reliability and accuracy of the results. A sufficient number of iterations is reached when the output no longer changes significantly with further increases.
#' @param seed An integer used to initialize the random number generator for reproducibility.
#'   Setting a seed ensures that the sequence of random numbers generated during the unmixing is reproducible. This is useful for debugging, testing, and comparing results across different runs.
#'
#' @return A list of data frames, where each data frame contains the predicted apportionments for a specific tracer. The last element of the list is a data frame containing the **Consistency Index (CI)** for each tracer.
#'
#' @details The function performs an individual tracer analysis to evaluate the conservativeness and discriminant capacity of each tracer. For each tracer, it constructs a determined system of linear equations by combining it with a minimal set of other tracers.
#'
#' There are two methods for completing this minimal set:
#' 1. The **"virtual" method** fabricates the remaining tracers by randomly generating values. This approach isolates the tracer of interest from the influence of other measured tracers.
#' 2. The **"random" method** randomly selects the remaining tracers from the available dataset, providing an assessment of how the tracer performs in combination with others.
#'
#' @references
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2020). Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. *Science of The Total Environment*, *720*, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537
#'
#' @export
individual_tracer_analysis <- function(data, completion_method = "virtual", iter = 5000, seed = 123456L)
{
  set.seed(seed)

  source_n <- nrow(inputSource(data))

  if (source_n == 2) ################################################################
  {
    sources  <- inputSource(data)
    mixtures <- inputMixture(data)
    tracer_n <- ncol(inputMixture(data))-1

    results <- list()
    for (i0 in 0:(tracer_n-1)) {
      if (completion_method == "virtual") {
        results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, iter, seed)
      } else if (completion_method == "random") {
        results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, iter, seed)
      } else {
        stop("Invalid value for 'completion_method' parameter. Must be 'virtual' or 'random'.")
      }
    }

    for (i1 in 1:length(results)) {
      results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
      results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
    }

    # Prepare data for the CI calculation
    ################################################################
    # Removing those out of 0 - 1
    mn_CI <- list() 
    for (i2 in 1:length(results)) {
      new_CI <- results[[i2]][ 
        results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
        results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1  , 
      ]

      mn_CI[[i2]] <-  colMeans(new_CI[, c(3:4)])
    }

    CI_new <- list()
    CI_new0 <- c()
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
      CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
      CI_new0[3]   <- 1 #tmp
      CI_new0[4]   <- 1 #tmp

      CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 4])) {
        CI_new[[i3]][, c(3, 4)] <- c(-1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }

    # New CI index calculation (distance equation)
    ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
        CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
        CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1  , ]
      prct_in[[i4]] <- nrow(tmp)/iter*100
    }
    tmp1 <- do.call("rbind", prct_in)

    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- 1 #tmp
      CI_new0[4]   <- 1 #tmp

      CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 4])) {
        CI_new[[i3]][, c(3, 4)] <- c(-1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else if (source_n == 3) ################################################################
  {
    sources  <- inputSource(data)
    mixtures <- inputMixture(data)
    tracer_n <- ncol(inputMixture(data))-1

    results <- list()
    for (i0 in 0:(tracer_n-1)) {
      if (completion_method == "virtual") {
        results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, iter, seed)
      } else if (completion_method == "random") {
        results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, iter, seed)
      } else {
        stop("Invalid value for 'completion_method' parameter. Must be 'virtual' or 'random'.")
      }
    }

    for (i1 in 1:length(results)) {
      results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
      results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
      results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
    }

    # Prepare data for the CI calculation
    ################################################################
    # Removing those out of 0 - 1
    mn_CI <- list() 
    for (i2 in 1:length(results)) {
      new_CI <- results[[i2]][ 
      results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
      results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
      results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1  , ]

      mn_CI[[i2]] <-  colMeans(new_CI[, c(3:5)])
    }

    CI_new <- list()
    CI_new0 <- c()
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
      CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
      CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
      CI_new0[4]   <- 1 #tmp
      CI_new0[5]   <- 1 #tmp

      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 5])) {
        CI_new[[i3]][, c(3, 4, 5)] <- c(-1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
 
    # New CI index calculation (distance equation)
    ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
        CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
        CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1 &
        CI_new[[i4]][,5] > 0 & CI_new[[i4]][,5] < 1  , ]
      prct_in[[i4]] <-  nrow(tmp)/iter*100
    }
    tmp1 <- do.call("rbind", prct_in)

    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- results[[i3]][5]
      CI_new0[4]   <- 1 #tmp
      CI_new0[5]   <- 1 #tmp

      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 5])) {
        CI_new[[i3]][, c(3, 4, 5)] <- c(-1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else if (source_n == 4) ################################################################
  {
    sources  <- inputSource(data)
    mixtures <- inputMixture(data)
    tracer_n <- ncol(inputMixture(data))-1

    results <- list()
    for (i0 in 0:(tracer_n-1)) {
      if (completion_method == "virtual") {
        results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, iter, seed)
      } else if (completion_method == "random") {
        results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, iter, seed)
      } else {
        stop("Invalid value for 'completion_method' parameter. Must be 'virtual' or 'random'.")
      }
    }

    for (i1 in 1:length(results)) {
      results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
      results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
      results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
      results[[i1]][, 6] <- as.numeric(as.character(results[[i1]][, 6]))
    }

    # Prepare data for the CI calculation
    ################################################################
    # Removing those out of 0 - 1
    mn_CI <- list() 
    for (i2 in 1:length(results)) {
      new_CI <- results[[i2]][ 
        results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
        results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
        results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1 &
        results[[i2]][,6] > -0.1 & results[[i2]][,6] < 1.1  , ]

      mn_CI[[i2]] <-  colMeans(new_CI[, c(3:6)])
    }

    CI_new <- list()
    CI_new0 <- c()
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
      CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
      CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
      CI_new0[4]   <- results[[i3]][6] - mn_CI[[i3]][4] + 1 / source_n
      CI_new0[5]   <- 1 #tmp
      CI_new0[6]   <- 1 #tmp
        
      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6])) {
          CI_new[[i3]][, c(3, 4, 5, 6)] <- c(-1, -1, -1, -1)
      } else {
          CI_new[[i3]] <- CI_new[[i3]]
      }
    }

    # New CI index calculation (distance equation)
    ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
        CI_new[[i4]][,3] > 0 & CI_new[[i4]][,3] < 1 &
        CI_new[[i4]][,4] > 0 & CI_new[[i4]][,4] < 1 &
        CI_new[[i4]][,5] > 0 & CI_new[[i4]][,5] < 1 &
        CI_new[[i4]][,6] > 0 & CI_new[[i4]][,6] < 1  , ]
      prct_in[[i4]] <-  nrow(tmp)/iter*100
    }
    tmp1 <- do.call("rbind", prct_in)

    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- results[[i3]][5]
      CI_new0[4]   <- results[[i3]][6]
      CI_new0[5]   <- 1 #tmp
      CI_new0[6]   <- 1 #tmp

      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6])) {
        CI_new[[i3]][, c(3, 4, 5, 6)] <- c(-1, -1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else if (source_n == 5) ################################################################
  {
    sources  <- inputSource(data)
    mixtures <- inputMixture(data)
    tracer_n <- ncol(inputMixture(data))-1

    results <- list()
    for (i0 in 0:(tracer_n-1)) {
      if (completion_method == "virtual") {
        results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, iter, seed)
      } else if (completion_method == "random") {
        results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, iter, seed)
      } else {
        stop("Invalid value for 'completion_method' parameter. Must be 'virtual' or 'random'.")
      }
    }

    for (i1 in 1:length(results)) {
      results[[i1]][, 3] <- as.numeric(as.character(results[[i1]][, 3]))
      results[[i1]][, 4] <- as.numeric(as.character(results[[i1]][, 4]))
      results[[i1]][, 5] <- as.numeric(as.character(results[[i1]][, 5]))
      results[[i1]][, 6] <- as.numeric(as.character(results[[i1]][, 6]))
      results[[i1]][, 7] <- as.numeric(as.character(results[[i1]][, 7]))
    }


    # Prepare data for the CI calculation
    ################################################################
    # Removing those out of 0 - 1
    mn_CI <- list() 
    for (i2 in 1:length(results)) {
      new_CI <- results[[i2]][ 
        results[[i2]][,3] > -0.1 & results[[i2]][,3] < 1.1 &
        results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1 &
        results[[i2]][,5] > -0.1 & results[[i2]][,5] < 1.1 &
        results[[i2]][,6] > -0.1 & results[[i2]][,6] < 1.1 &
        results[[i2]][,7] > -0.1 & results[[i2]][,7] < 1   , ]
      mn_CI[[i2]] <-  colMeans(new_CI[, c(3:7)])
    }

    CI_new <- list()
    CI_new0 <- c()
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
      CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
      CI_new0[3]   <- results[[i3]][5] - mn_CI[[i3]][3] + 1 / source_n
      CI_new0[4]   <- results[[i3]][6] - mn_CI[[i3]][4] + 1 / source_n
      CI_new0[5]   <- results[[i3]][7] - mn_CI[[i3]][5] + 1 / source_n
      CI_new0[6]   <- 1 #tmp
      CI_new0[7]   <- 1 #tmp
      CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)]

      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6]|is.na(CI_new[[i3]][1, 7]))) {
        CI_new[[i3]][, c(3, 4, 5, 6, 7)] <- c(-1, -1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }

    # New CI index calculation (distance equation)
    ################################################################
    prct_in <- list()
    for (i4 in 1:length(results)) {
      tmp <- CI_new[[i4]][ 
         CI_new[[i4]][, 3] > 0 & CI_new[[i4]][, 3] < 1 &
         CI_new[[i4]][, 4] > 0 & CI_new[[i4]][, 4] < 1 &
         CI_new[[i4]][, 5] > 0 & CI_new[[i4]][, 5] < 1 &
         CI_new[[i4]][, 6] > 0 & CI_new[[i4]][, 6] < 1 &
         CI_new[[i4]][, 7] > 0 & CI_new[[i4]][, 7] < 1  , ]

      prct_in[[i4]] <-  nrow(tmp) / iter * 100
    }
    tmp1 <- do.call("rbind", prct_in)
    
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- results[[i3]][5]
      CI_new0[4]   <- results[[i3]][6]
      CI_new0[5]   <- results[[i3]][7]
      CI_new0[6]   <- 1 #tmp
      CI_new0[7]   <- 1 #tmp
      CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)]
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6]|is.na(CI_new[[i3]][1, 7]))) {
        CI_new[[i3]][, c(3, 4, 5, 6, 7)] <- c(-1, -1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else {
    stop(paste0("Error: individual_tracer_analysis is not implemented for ", source_n," sources."))
  }

  # Prepare for output
  ################################################################
  names_CI_R <- colnames(inputMixture(data))[-c(1)]
  names_CI_R <- as.data.frame(names_CI_R)
  CI <- cbind(names_CI_R, tmp1)
  CI_new$DB <- data # add the database at the end
  CI_new$CI <- CI   # add the CI index at the end

  return(CI_new)
}

