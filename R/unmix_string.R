#' Unmix sediment mixtures
#'
#' Asses the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#'
#' @param data Data frame containing sediment source and mixtures
#' @param samples Number of samples in each hypercube dimension
#' @param iter Iterations in the source variability analysis
#' @param Means Boolean to switch when using mean and sd data
#' @param seed Seed for the random number generator
#'
#' @return Data frame containing the relative contribution of the sediment sources for each sediment mixture and iterations
#'
#' @export
#'
unmix_string <- function(data, samples = 100L, iter = 100L, Means = F, seed = 123456L){
  system.time({
    if (Means == T) {
      sources <- data[c(1:(nrow(data)-1)),c(2:ncol(data))]
      x <- round((ncol(data)-3)/2 + 2)
      mixtures <- inputSample(data[nrow(data),1:x])
      colnames(sources)[1] <- "id"
      
    } else {
      ###########################################
      sources <- inputSource(data)
      mixtures <- inputSample(data)
      ###########################################
    }
    
    # verify the number of sources and properties
    if (nrow(sources) -1 >= ncol(mixtures) ) {
      warning("As a minimum, n - 1 properties are required to discriminate rigorously between n sources. Additional properties are frequently required to increase the reliability of the results")
    }
    
    nsources1<- as.data.frame(unique(data[,2]))
    nsources<- nsources1[-nrow(nsources1),] 
    nsources<- as.vector(nsources) 
    
    if (Means == T) {
      cat("Summary of the model imputs:
        ", (ncol(data)-3)/2, "variables from",nrow(nsources1)-1,"sources (",nsources,")",
          "\n")

    } else {
      cat("Summary of the model imputs:
        ", ncol(data)-2, "variables from",nrow(nsources1)-1,"sources (",nsources,")",
          "\n")
    }
    
    #invisible(readline(prompt="Press [enter] to unmix your data"))
    
    # results <- unmix_c(sources, mixtures, samples, 1, seed)  
    if (iter==1) {  
      results <- unmix_c(sources, mixtures, samples, iter, seed)
    }
    else {  
      results <- unmix_c(sources, mixtures, samples, iter+1, seed)
    }
    
    if (nrow(sources) == 3) {
        results <- results[results[,3] >= 0.0, ]
        results <- results[results[,3] <= 1.0, ]
        results <- results[results[,4] >= 0.0, ]
        results <- results[results[,4] <= 1.0, ]
        results <- results[results[,5] >= 0.0, ]
        results <- results[results[,5] <= 1.0, ]
    } else if (nrow(sources) == 4) {
        results <- results[results[,3] >= 0.0, ]
        results <- results[results[,3] <= 1.0, ]
        results <- results[results[,4] >= 0.0, ]
        results <- results[results[,4] <= 1.0, ]
        results <- results[results[,5] >= 0.0, ]
        results <- results[results[,5] <= 1.0, ]
        results <- results[results[,6] >= 0.0, ]
        results <- results[results[,6] <= 1.0, ]
    } else if (nrow(sources) == 5) {
        results <- results[results[,3] >= 0.0, ]
        results <- results[results[,3] <= 1.0, ]
        results <- results[results[,4] >= 0.0, ]
        results <- results[results[,4] <= 1.0, ]
        results <- results[results[,5] >= 0.0, ]
        results <- results[results[,5] <= 1.0, ]
        results <- results[results[,6] >= 0.0, ]
        results <- results[results[,6] <= 1.0, ]
        results <- results[results[,7] >= 0.0, ]
        results <- results[results[,7] <= 1.0, ]
    } else {
      results <- results
    }

    
    cat('\n')
    cat(crayon::bold("Attention", "from a total of", iter, " only ", nrow(results), "solutions remain inside the physical space"))
    cat('\n')
    cat('\n')
    
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
    
    
    {
      if (iter==1) {  
        cat("Summary of the model outputs:",
            "\n",
            "See below the result/s of the unmixing process using the central value or the average with no correction",
            "\n",
            "\n")
        print(aggregate(. ~ id, data = results, function(x) c(mean = mean(x))))
      }  
      
      else {  
        cat("Summary of the model outputs:",
            "\n",
            "See below the result/s of the unmixing process using the source variability of the best", iter, "results, notice that the first row of the results is the central value or the average with no correction",
            "\n",
            "\n")
        print(aggregate(. ~ id, data = results, function(x) c(mean = mean(x), SD = sd(x))))
        
      }  
    }
    
    return(results)

  })
}
