#' CR values of each sediment mixture for four sources
#'
#' Compute the consensus ranking method for the selected dataset.
#'
#' @param source Data frame containing the sediment sources from a dataset
#' @param mixture Data frame containing one of the dataset mixtures
#' @param maxiter Number of iteration for each tracer
#' @param seed Seed for the random number generator
#' 
#' @return Data frame containing the CR score for each tracer.
#'
#' @export
#'
cr_ns <- function(source, mixture, maxiter = 2000, seed = 123456)
{
  set.seed(seed)
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
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
  pb <- txtProgressBar(min = 0, max = maxiter, style = 3, width = 50, char = "=")
  
  while (min(cr1) < maxiter)
  {
    var <- sample(c(1:cols), nsource+1)
    
    if(min(cr1[var]) < maxiter)
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
        if(cr1[var[i]] < maxiter)
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
            if(cr1[var[i]] < maxiter)
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
    setTxtProgressBar(pb, min(cr1))}
  
  score <- 100-((cr2*100)/cr1)
  cr <- data.frame(tracer, score)
  cr <- cr[rev(order(cr$score)),]
  row.names(cr) <- NULL
  cat('\n')
  cat('\n')
  return(cr)
}

