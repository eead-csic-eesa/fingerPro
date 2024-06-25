#' Conservativeness Index (CI)
#'
#' The function quantify the predictions of each individual tracer in your dataset and calculate the CI value of each tracer for 2S, 3S and 4S.
#'
#' @param data Data frame containing sediment source and mixtures
#' @param points Number of solutions per tracer
#' @param Means Boolean to switch when using mean and sd data
#' @param seed Seed for the random number generator
#'
#' @return List of data frames containing the selected number of possible prediction for each tracer. The last lists correspond to the dataset and the CI index.
#'
#' @export
#'  
CI_Method <- function(data, points = 2000, Means = T, triangles = "virtual", seed = 123456L)
 {
  set.seed(seed)
  
  source_n <- nrow(inputSource(data))
  system.time({
    
    if (source_n == 2) ###########################################################################################
    {
      if (Means == T) {
        sources1 <- data[-nrow(data), c(3:ncol(data))]
        id_sources <- (1:(nrow(data) - 1))
        id <- paste('S', id_sources, sep = '')
        sources <- cbind(id, sources1)
        mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
        
        results <- list()
        for (i0 in 0:(((ncol(data) + 1) / 2) - 3)) {
            if (triangles == "virtual") {
              results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
            } else if (triangles == "random") {
              results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
            } else {
              warning("Not applicable")
            }
        }
      } else {
        ############################################
        sources  <- inputSource(data)
        mixtures <- inputSample(data)
        ###########################################
        results <- list()
        for (i0 in 0:(ncol(data) - 3)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
          } else if (triangles == "random") {
            results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
          } else {
            warning("Not applicable")
          }
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
            results[[i2]][,4] > -0.1 & results[[i2]][,4] < 1.1  , ]

        mn_CI[[i2]] <-  colMeans(new_CI[, c(3:4)])
      }
      
      CI_new <- list()
      CI_new0 <- c()
      for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3] - mn_CI[[i3]][1] + 1 / source_n
        CI_new0[2]   <- results[[i3]][4] - mn_CI[[i3]][2] + 1 / source_n
        CI_new0[3]   <- 1 #tmp
        CI_new0[4]   <- 1 #tmp
        
        CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)] # forma cutre
        
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
        prct_in[[i4]] <-  nrow(tmp)/points*100
      }
      tmp1 <- do.call("rbind", prct_in)
      
      for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3]
        CI_new0[2]   <- results[[i3]][4]
        CI_new0[3]   <- 1 #tmp
        CI_new0[4]   <- 1 #tmp
        
        CI_new[[i3]] <- CI_new0[, c(3, 4, 1, 2)] # forma cutre
        
        # Loop to prevent NA when no solutions are inside the range
        if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 4])) {
          CI_new[[i3]][, c(3, 4)] <- c(-1, -1)
        } else {
          CI_new[[i3]] <- CI_new[[i3]]
        }
      }
      
    } else if (source_n == 3) ###########################################################################################
    {
      if (Means == T) {
        sources1 <- data[-nrow(data), c(3:ncol(data))]
        id_sources <- (1:(nrow(data) - 1))
        id <- paste('S', id_sources, sep = '')
        sources <- cbind(id, sources1)
        mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
        
        results <- list()
        for (i0 in 0:(((ncol(data) + 1) / 2) - 3)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
          } else if (triangles == "random") {
            results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
          } else {
            warning("Not applicable")
          }
        }
      } else {
        ############################################
        sources  <- inputSource(data)
        mixtures <- inputSample(data)
        ###########################################
        results <- list()
        for (i0 in 0:(ncol(data) - 3)) {
          if (triangles == "virtual") {
            results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
          } else if (triangles == "random") {
            results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
          } else {
            warning("Not applicable")
          }
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
      
      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)] # forma cutre
      
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
      prct_in[[i4]] <-  nrow(tmp)/points*100
    }
    tmp1 <- do.call("rbind", prct_in)
    
    for (i3 in 1:length(results)) {
      CI_new0      <- results[[i3]][3]
      CI_new0[2]   <- results[[i3]][4]
      CI_new0[3]   <- results[[i3]][5]
      CI_new0[4]   <- 1 #tmp
      CI_new0[5]   <- 1 #tmp
      
      CI_new[[i3]] <- CI_new0[, c(4, 5, 1, 2, 3)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 5])) {
        CI_new[[i3]][, c(3, 4, 5)] <- c(-1, -1, -1)
      } else {
        CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
  } else if (source_n == 4)	###########################################################################################
    {
    if (Means == T) {
      sources1 <- data[-nrow(data), c(3:ncol(data))]
      id_sources <- (1:(nrow(data) - 1))
      id <- paste('S', id_sources, sep = '')
      sources <- cbind(id, sources1)
      mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
      
      results <- list()
      for (i0 in 0:(((ncol(data) + 1) / 2) - 3)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
        } else if (triangles == "random") {
          results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
        } else {
          warning("Not applicable")
        }
      }
    } else {
      ############################################
      sources  <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
      results <- list()
      for (i0 in 0:(ncol(data)-3)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
        } else if (triangles == "random") {
          results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
        } else {
          warning("Not applicable")
        }
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
        
      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)] # forma cutre
      
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
      prct_in[[i4]] <-  nrow(tmp)/points*100
    }
    tmp1 <- do.call("rbind", prct_in)

    for (i3 in 1:length(results)) {
        CI_new0      <- results[[i3]][3]
        CI_new0[2]   <- results[[i3]][4]
        CI_new0[3]   <- results[[i3]][5]
        CI_new0[4]   <- results[[i3]][6]
        CI_new0[5]   <- 1 #tmp
        CI_new0[6]   <- 1 #tmp
        
      CI_new[[i3]] <- CI_new0[, c(5, 6, 1, 2, 3, 4)] # forma cutre
      
      # Loop to prevent NA when no solutions are inside the range
      if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6])) {
          CI_new[[i3]][, c(3, 4, 5, 6)] <- c(-1, -1, -1, -1)
      } else {
          CI_new[[i3]] <- CI_new[[i3]]
      }
    }
    
    } else if (source_n == 5)	####################################################################################################
      {
    if (Means == T) {
      sources1 <- data[-nrow(data), c(3:ncol(data))]
      id_sources <- (1:(nrow(data) - 1))
      id <- paste('S', id_sources, sep = '')
      sources <- cbind(id, sources1)
      mixtures <- inputSample(data[, c(1:((ncol(data) / 2) + 1))])
      
      results <- list()
      for (i0 in 0:(((ncol(data) + 1) / 2) - 3)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
        } else if (triangles == "random") {
          results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
        } else {
          warning("Not applicable")
        }
      }
    } else {
      ############################################
      sources <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
      results <- list()
      for (i0 in 0:(ncol(data) - 3)) {
        if (triangles == "virtual") {
          results[[i0 + 1]] <- triangles_virtual_c(sources, mixtures, i0, points, seed)
        } else if (triangles == "random") {
          results[[i0 + 1]] <- triangles_random_c(sources, mixtures, i0, points, seed)
        } else {
          warning("Not applicable")
        }
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
    CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)] # forma cutre
    
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
  
      prct_in[[i4]] <-  nrow(tmp) / points * 100
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
    CI_new[[i3]] <- CI_new0[, c(6, 7, 1, 2, 3, 4, 5)] # forma cutre
    
    # Loop to prevent NA when no solutions are inside the range
    if (is.na(CI_new[[i3]][1, 3])|is.na(CI_new[[i3]][1, 6]|is.na(CI_new[[i3]][1, 7]))) {
      CI_new[[i3]][, c(3, 4, 5, 6, 7)] <- c(-1, -1, -1, -1)
    } else {
      CI_new[[i3]] <- CI_new[[i3]]
    }
  }
    
   } else 
     {
        print("ERROR: not implemented, please, contact code developers || lizaga.ivan10@gmail.com")
      }
  
  # Prepare for output
  ################################################################
      if (Means == T) {
        names_CI_R <- names(data[3:ceiling((ncol(data)/2))])# round function do 0.5 = 0
      }
      if (Means == F) {
        names_CI_R <- names(data[3:ncol(data)])
      }
      
      names_CI_R <- as.data.frame(names_CI_R)
      CI <- cbind(names_CI_R, tmp1)
      CI_new$DB <- data # add the database at the end
      CI_new$CI <- CI   # add the CI index at the end
      
      return(CI_new)
      
  })
}

