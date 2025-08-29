#' Extract all possible tracer combinations with one tracer.
#'
#' This function generates a list of all possible tracer combinations to identify the most discriminant and serves as a seed to build a consistent tracer selection in a subsequent step. This analysis explores minimal tracer combinations (one tracer for two sources) and solves the resulting determined system of equations to assess the variability of each combination. The dispersion of the solution reflects the discriminant capacity of each tracer combination: a lower dispersion indicates a higher discriminant capacity. Typically, the most discriminant tracer combination corresponds to the result of DFA analysis. In this analysis, the solutions are not restricted to the physically feasible space, which can be valuable for identifying problematic tracer selections that might be masked when using constrained unmixing models.
#'
#' @param source Data frame containing the sediment sources from a dataset.
#' @param mixture Data frame containing one of the dataset mixtures.
#' @param iter Iterations in the variability analysis of each tracer combination.
#' @param seed An integer value used to initialize the random number generator.
#'   Setting a seed ensures that the sequence of random numbers generated during the unmixing is reproducible. This is useful for debugging, testing, and comparing results across different runs.
#'   If no seed is provided, a random seed will be generated.
#' 
#' @return A data frame containing all possible tracer combinations from the dataset. Each combination is characterized by its corresponding average solution and dispersion (standard deviation), as well as the percentage of solutions that fall within the physically feasible space.
#'
CTS_seeds_singles <- function(source, mixture, iter = 1000, seed = 123456)
{
  set.seed(seed)
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
  cols <- (ncol(source)-1)/2
  tracer <- colnames(source)[1:cols]
  
  n <- cols*2+1 # n column
  
  df <- data.frame(tracers=character(), w1=double(), w2=double(), sd_w1=double(), sd_w2=double(), percent_physical=double(), stringsAsFactors=FALSE)
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = cols, style = 3, width = 50, char = "=")
  
  for (i in c(1:cols))
  {
    w1 <- c()
    w2 <- c()
    sc <- 0
    snc <- 0
    
    for (l in c(1:iter))
    {
      if(l == 1)
      {
        s11 <- source[1, i]
        s21 <- source[2, i]
        m1  <- mixture[i]
      }
      else
      {
        s11 <- source[1, i] + source[1, i+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
        s21 <- source[2, i] + source[2, i+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
        m1  <- mixture[i]
      }
      
      # maxima script
      #matrix(
      # [1,1],
      # [s11,s21]
      #);
      #M:%;
      #determinant(M);
      #invert(M)*determinant(M);
      #ratsimp(%);
      
      det = s21-s11
      
      if(det != 0.0)
      {
        x1 =  s21 - m1
        x2 = -s11 + m1
        
        x1 = x1 / det
        x2 = x2 / det

        if(min(x1,x2) >= 0.0 && max(x1,x2) <= 1.0)
        {
          sc <- sc + 1
        }
        else
        {
          snc <- snc + 1
        }
        
        w1 <- c(w1, x1)
        w2 <- c(w2, x2)
      }
    }
    df[nrow(df) + 1,] <- list(paste0(tracer[i]), w1[1], w2[1], (quantile(w1, 0.84)[[1]] - quantile(w1, 0.16)[[1]])/2, (quantile(w2, 0.84)[[1]] - quantile(w2, 0.16)[[1]])/2, sc/(sc+snc))
		setTxtProgressBar(pb, i)
  }
  
  df$max_sd_wi <- pmax(df$sd_w1, df$sd_w2)
  df <- df[order(df$max_sd_wi),]
  row.names(df) <- NULL
  cat('\n')
  cat('\n')
  return(df)
}

