#' @title Consensus Ranking (CR) method for tracer selection
#'
#' @description This function computes the Consensus Ranking (CR) method, an ensemble technique to identify non-conservative and dissenting tracers in sediment fingerprinting studies. The method combines predictions from single-tracer models and is based on a scoring function derived from a series of random "debates" between tracers.
#'
#' @param data A data frame containing sediment source and mixture data. Users should ensure their data is in a valid format by using the `check_database()` function before running the CR method.
#' @param debates An integer specifying the target number of debates each tracer should participate in. The function will run until each tracer has participated in at least this many debates.
#' @param seed An integer used to initialize the random number generator for reproducibility.
#'
#' @return A data frame containing the CR score for each tracer. The score, ranging from 100 to 0, indicates the tracer's rank in terms of consensus and conservativeness. Tracers are ordered by their score in descending order, with the most conservative tracers having high scores and dissenting tracers having low scores.
#'
#' @details The Consensus Ranking method is based on a series of random debates to test the compatibility of tracers. In each debate, a random subset of tracers is selected. The size of this subset is determined by the number of sources, corresponding to the minimum number of equations needed to overdetermine the unmixing model.
#'
#' For each debate, a least-squares method is used to find a solution to the overdetermined mass balance equations. The consensus of the debate is measured by the mathematical compatibility of the tracers, specifically using the Root Mean Square Error (RMSE) of the mass balance equations. The tracer whose exclusion from the debate results in lowest RMSE is identified as the "dissenting" tracer for that round.
#'
#' This process is repeated for a specified number of debates. Each tracer accumulates a count of total participations and a count of lost debates (being identified as dissenting). The final CR score is a quantitative measure of consensus, calculated as `100 - (lost debates / total debates) * 100`.
#'
#' A low CR score indicates that a tracer frequently disrupts the consensus and is considered a non-conservative or dissenting tracer. Conversely, a high CR score suggests the tracer is in frequent agreement with the others, making it a reliable and conservative tracer for the unmixing model. This method is robust and does not require pre-screening or filtering of tracers.
#'
#' @references
#' Lizaga, I., Latorre, B., Gaspar, L., & Navas, A. (2020). Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. *Science of The Total Environment*, *720*, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537
#'
#' @export
CR <- function(data, debates = 1000, seed = 123456)
{
	source <- inputSource(data)
	mixture <- inputMixture(data)
	
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])

  set.seed(seed)
    
  # normalize
  cols <- (ncol(source)-1)/2
  for (col in c(1:cols))
  {
    mx <- max(source[,col]+source[,cols+col])
    mn <- min(source[,col]-source[,cols+col])
    source[,col] <- ( source[,col] - mn ) / ( mx - mn )
    source[,cols+col] <- source[,cols+col] / ( mx - mn )
    mixture[,col] <- ( mixture[,col] - mn ) / ( mx - mn )
  }
  
  tracer <- colnames(source)[1:cols]
  tracer <- gsub("^mean_", "", tracer)
  
  cr1 <- c()
  cr2 <- c()
  for (col in c(1:cols))
  {
    cr1[col] <- 0
    cr2[col] <- 0
  }
  
  nsource <- nrow(source)
  ntracer <- (ncol(source)-1)/2
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = debates, style = 3, width = 50, char = "=")
  
  while (min(cr1) < debates)
  {
    var <- sample(c(1:cols), nsource+1)
    
    if(min(cr1[var]) < debates)
    {
      gof <- c()
      max <- c()
      min <- c()
      rnd <- matrix(0, nrow = nsource, ncol = nsource+1)
      for (i in c(1:(nsource+1)))
      {
        j <- cols*2+1 # n column
         for (k in c(1:nsource))
         {
      	  rnd[k, i] <- rt(1, source[k, j]) / sqrt(source[k, j])
      	 }
      }
      for (i in c(1:(nsource+1)))
      {		
        j <- c(1:(nsource+1))[-i]

      	y <- c()
      	x <- matrix(, nrow = nsource, ncol = nsource-1)
    		csource <- matrix(, nrow = nsource, ncol = nsource)
          	
        ls <- c()
        for (jj in c(1:nsource))
        {
    	    k <- j[jj]
     	   	l <- var[k]
          csource[nsource, jj] <- source[nsource, l] + source[nsource, ntracer+l] * rnd[nsource, k]
          ls <- c(ls, csource[nsource, jj])

          y <- c(y, mixture[l] - ls[jj])
          for (kk in c(1:(nsource-1)))
          {
            csource[kk,jj] <- source[kk,l] + source[kk, ntracer+l] * rnd[kk, k]
            x[jj, kk] <- csource[kk,jj] - ls[jj]
          }
        }

		    # least squares method
		    model <- lm.fit(x, y)
		    x1 <- as.vector(coef(model))
		    x1 <- c(x1, 1-sum(x1))

				# compute gof
				err <- 0
        for (jj in c(1:nsource))
        {
    	    k <- j[jj]
     	   	l <- var[k]
        	vm <- 0
        	for (kk in c(1:nsource))
          {
          		vm <- vm + csource[kk,jj] * x1[kk]
          }
          err <- err + (vm - mixture[l])*(vm - mixture[l])
        }

        gof <- c(gof, err)
        max <- c(max, max(x1))
        min <- c(min, min(x1))
      }
      
      sgof <- sort(gof)
      if(sgof[nsource+1]-sgof[nsource] > 10.0*(sgof[nsource]-sgof[1]))
      {
        i <- which.max(gof)
        if(cr1[var[i]] < debates)
        {
          cr1[var[i]] <- cr1[var[i]] + 1
          cr2[var[i]] <- cr2[var[i]] + 1
        }
      }
      else
      {
        win <- which.min(gof)
        if(max[win] <= 1.0 && min[win] >= 0.0)
        {
          for (i in c(1:(nsource+1)))
          {
            if(cr1[var[i]] < debates)
            {
              cr1[var[i]] <- cr1[var[i]] + 1
              
              if(i == win)
              {
                cr2[var[i]] <- cr2[var[i]] + 1
              }
            }
          }
        }
      }
    }
    setTxtProgressBar(pb, min(cr1))
  }
  
  CR_score <- 100-((cr2*100)/cr1)
  cr <- data.frame(tracer, CR_score)
  cr <- cr[rev(order(cr$CR_score)),]
  row.names(cr) <- NULL
  cat('\n')
  cat('\n')
  return(cr)
}

