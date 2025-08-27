#' Singles function 
#'
#' Extract all the possible selections formed by sets of n ??? 1 tracers and solve them by using standard methods
#'
#' @param source Data frame containing the sediment sources from a dataset
#' @param mixture Data frame containing one of the dataset mixtures
#' @param iter Number of iteration for each tracer
#' @param seed Seed for the random number generator
#' 
#' @return Data frame containing all the possible combination from your dataset, their system of equation solution, the consistency and the discriminant capacity. 
#'
#' @export
#'
singles <- function(source, mixture, iter = 1000, seed = 123456)
{
  {
  set.seed(seed)
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
  cols <- (ncol(source)-1)/2
  tracer <- colnames(source)[1:cols]
  
  n <- cols*2+1 # n column
  
  df <- data.frame(id=character(), w1=double(), w2=double(), Dw1=double(), Dw2=double(), cons=double(), stringsAsFactors=FALSE)
  
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
  }
  setTxtProgressBar(pb, i)}
  
  df <- transform(df, Dmax = pmax(Dw1, Dw2))
  df <- df[order(df$Dmax),]
  row.names(df) <- NULL
  cat('\n')
  cat('\n')
  return(df)
}
