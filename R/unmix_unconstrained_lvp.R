#' Unmix sediment mixtures using unconstrained (least-squares) optimization and linear variability propagation (LVP)
#'
#' Asses the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#'
#' @param data Data frame containing sediment source and mixtures
#' @param variability Character string specifying the type of variability to calculate.
#'   Possible values are "SD" for Standard Deviation or "SEM" for Standard Error of the Mean.
#' @param iter Iterations in the source variability analysis.
#' @param means Boolean to switch when using mean and sd data
#' @param seed Seed for the random number generator
#'
#' @return Data frame containing the relative contribution of the sediment sources for each sediment mixture and iterations
#'
#' @export
#'

unmix_unconstrained_lvp <- function(data, variability = "SEM", iter = 1000, means = F, seed = 123456L)
{
  set.seed(seed)
  
  sr_names <- unique(data[, 2])

#  if(is.numeric(iter)){ iter <- iter} else { 
#    if(iter == "very short") iter <- 500
#    if(iter == "short")      iter <- 1000
#    if(iter == "medium")     iter <- 5000
#    if(iter == "long")       iter <- 20000
#    if(iter == "bananas")    iter <- 100000
#  }
  
  if (means == T) {
    source <- data[c(1:(nrow(data) - 1)), c(2:ncol(data))]
    x <- round((ncol(data) - 3) / 2 + 2)
    mixture <- inputMixture(data[nrow(data), 1:x])
    colnames(source)[1] <- "id"
    colnames(mixture) <- colnames(source[, c(1:(ncol(source) / 2))])
    
  } else {
    ###########################################
    source <- inputSource(data)
    mixture <- inputMixture(data)
    ###########################################
  }
  
  if (nrow(source)-1 >= ncol(mixture)) { #ntracer >= nsources - 1
    warning("As a minimum, n - 1 properties are required to discriminate rigorously between n sources. Additional properties are frequently required to increase the reliability of the results")
  }

    snames <- source[,1]
    
    source <- data.matrix(source[-1])
    mixture <- data.matrix(mixture[-1])
    
    # normalize
    cols <- (ncol(source)-1)/2
    for (col in c(1:cols))
    {
      mx <- max(source[, col] + source[, cols + col])
      mn <- min(source[, col] - source[, cols + col])
      source[, col] <- (source[, col] - mn) / (mx - mn)
      source[, cols + col] <- source[, cols + col] / (mx - mn)
      mixture[, col] <- (mixture[, col] - mn) / (mx - mn)
    }
    
    nsource <- nrow(source)
    ntracer <- (ncol(source)-1)/2
    
    r <- matrix(, nrow = iter, ncol = nsource+2)
    
    # compute central solution
    y <- c()
    x <- matrix(, nrow = ntracer, ncol = nsource-1)
    for (j in c(1:ntracer))
    {
      y <- c(y, mixture[1,j][[1]]-source[nsource,j])
      for (k in c(1:(nsource-1)))
      {
        x[j, k] <- source[k,j]-source[nsource,j]
      }
    }
    
    cat(crayon::cyan(crayon::bold("Summary of the model inputs:\n")))
    cat(ntracer, "variables from", nsource, "sources", "\n")
    
    
    # Introduce the progress bar
    pb <- txtProgressBar(min = 0, max = iter, style = 3, width = 50, char = "=")
    
    
    # least squares method
    model <- lm.fit(x, y)
    cw <- as.vector(coef(model))
    cw <- c(cw, 1-sum(cw))
    
    for (i in c(1:iter))
    {
      vm <- c(1:ntracer)*0		
      if(i==2 || i==3)
      {
        # use measured mixture
        for (j in c(1:ntracer))
        {
          vm[j] <- mixture[1,j][[1]]
        }
      }
      else
      {
        # compute virtual mixture
        for (j in c(1:ntracer))
        {
          vm[j] <- 0
          for (k in c(1:nsource))
          {
            
            if (variability == "SEM") {
              x1 <- rt(1, source[k, ntracer*2+1]) / sqrt(source[k, ntracer*2+1])
              
            } else if (variability == "SD") {
              x1 <- rt(1, source[k, ntracer*2+1])
              
            } else  {
              cat("Not implemented. Please contact code developers")
              
            }
            
            vm[j] <- vm[j] + cw[k] * (source[k,j] + source[k,ntracer+j] * x1 )
            
           }
        }
      }
      
      y <- c()
      x <- matrix(, nrow = ntracer, ncol = nsource-1)
      for (j in c(1:ntracer))
      {
        y <- c(y, vm[j]-source[nsource,j])
        for (k in c(1:(nsource-1)))
        {
          x[j, k] <- source[k,j]-source[nsource,j]
        }
      }
      # least squares method
      model <- lm.fit(x, y)
      w <- as.vector(coef(model))
      w <- c(w, 1-sum(w))
      
      gof <- c()
      for (j in c(1:ntracer))
      {
        x1 <- 0.0
        for (k in c(1:nsource))
        {
          x1 <- x1 + w[k] * source[k,j]
        }
        # gof <- c(gof, (vm[j]-x1)^2)
        gof <- c(gof, (mixture-x1)^2)
      }
      
      gof <- 1.0 - mean(gof)
      gof <- max(gof, 0)
      
      w <- c(1, gof, w)
      
      r[i, ] <- w
      
      setTxtProgressBar(pb, min(i))
    }
    
    r <- as.data.frame(r)
    # colnames(r) <- c('id', 'gof', paste( "w.", snames[c(1:nsource)], sep=""))
    colnames(r) <- c('id', 'gof', paste( sr_names[c(1:nsource)], sep=""))
    
    r2 <- r
    r2 <- r2[, c(setdiff(colnames(r2), "gof"))]
    
    cat("\n", "\n")
    cat(crayon::green(crayon::bold("Summary of the model outputs:")))
    cat("\n", 
        "The results below show the outcome of the unmixing process, which incorporated a source variability analysis \nover", 
        iter, "iterations. Notice that the second and third rows of the result correspond to the solution for the \ncentral or mean value of the sources.", 
        "\n", "\n")
    
     print(aggregate(. ~ id, data = r2, function(x) c(mean = mean(x), SD = sd(x))))
     
     
     invisible(return(r))
  }
