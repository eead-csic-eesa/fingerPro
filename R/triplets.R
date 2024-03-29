#' triplets function 
#'
#' Extract all the possible selections formed by sets of n ??? 1 tracers and solve them by using standard methods
#'
#' @param source Data frame containing the sediment sources from a dataset
#' @param mixture Data frame containing one of the dataset mixtures
#' @param iter Number of iteration for each tracer
#' @param seed Seed for the random number generator
#' 
#' @return Data frame containing all the possible triplets combination from your dataset, their system of equation solution, the consistency and the discriminant capacity. 
#'
#' @export
#'
triplets <- function(source, mixture, iter = 1000, seed = 123456)
{
  set.seed(seed)
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
  cols <- (ncol(source)-1)/2
  tracer <- colnames(source)[1:cols]
  
  n <- cols*2+1 # n column
  
  df <- data.frame(id=character(), w1=double(), w2=double(), w3=double(), w4=double(), Dw1=double(), Dw2=double(), Dw3=double(), Dw4=double(), cons=double())
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = cols, style = 3, width = 50, char = "=")
  
  for (i in c(1:cols))
  {
    for (j in c(i+1:cols))
    {
      for (k in c(j+1:cols))
      {
        if(i<j & j<k & k<=cols)
        {
          w1 <- c()
          w2 <- c()
          w3 <- c()
          w4 <- c()
          sc <- 0
          snc <- 0
          
          for (l in c(1:iter))
          {
            if(l == 1)
            {
              s11 <- source[1, i]
              s21 <- source[2, i]
              s31 <- source[3, i]
              s41 <- source[4, i]
              m1  <- mixture[i]
              
              s12 <- source[1, j]
              s22 <- source[2, j]
              s32 <- source[3, j]
              s42 <- source[4, j]
              m2  <- mixture[j]
              
              s13 <- source[1, k]
              s23 <- source[2, k]
              s33 <- source[3, k]
              s43 <- source[4, k]
              m3  <- mixture[k]
            }
            else
            {
              s11 <- source[1, i] + source[1, i+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
              s21 <- source[2, i] + source[2, i+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
              s31 <- source[3, i] + source[3, i+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
              s41 <- source[4, i] + source[4, i+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
              m1  <- mixture[i]
              
              s12 <- source[1, j] + source[1, j+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
              s22 <- source[2, j] + source[2, j+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
              s32 <- source[3, j] + source[3, j+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
              s42 <- source[4, j] + source[4, j+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
              m2  <- mixture[j]
              
              s13 <- source[1, k] + source[1, k+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
              s23 <- source[2, k] + source[2, k+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
              s33 <- source[3, k] + source[3, k+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
              s43 <- source[4, k] + source[4, k+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
              m3  <- mixture[k]
            }
            
            det <- s21*(s32*s43-s33*s42)-s11*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+s11*(s22*s43-s23*s42)+s31*(s12*s43-s13*s42)-s21*(s12*s43-s13*s42)+(s22*s33-s23*s32)*s41-(s12*s33-s13*s32)*s41+
              (s12*s23-s13*s22)*s41-s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31
            
            if(det != 0.0)
            {
              x1 = s21*(s32*s43-s33*s42)-s31*(s22*s43-s23*s42)+(s22*s33-s23*s32)*s41 + m1*(-s32*s43+s22*s43+s33*s42-s23*s42-s22*s33+s23*s32) + m2*(s31*s43-s21*s43-s33*s41+s23*s41+s21*s33-s23*s31) + m3*(-s31*s42+s21*s42+s32*s41-s22*s41-s21*s32+s22*s31)
              
              x2 = -s11*(s32*s43-s33*s42)+s31*(s12*s43-s13*s42)-(s12*s33-s13*s32)*s41 + m1*(s32*s43-s12*s43-s33*s42+s13*s42+s12*s33-s13*s32) + m2*(-s31*s43+s11*s43+s33*s41-s13*s41-s11*s33+s13*s31) + m3*(s31*s42-s11*s42-s32*s41+s12*s41+s11*s32-s12*s31)
              
              x3 = s11*(s22*s43-s23*s42)-s21*(s12*s43-s13*s42)+(s12*s23-s13*s22)*s41 + m1*(-s22*s43+s12*s43+s23*s42-s13*s42-s12*s23+s13*s22) + m2*(s21*s43-s11*s43-s23*s41+s13*s41+s11*s23-s13*s21) + m3*(-s21*s42+s11*s42+s22*s41-s12*s41-s11*s22+s12*s21)
              
              x4 = -s11*(s22*s33-s23*s32)+s21*(s12*s33-s13*s32)-(s12*s23-s13*s22)*s31 + m1*(s22*s33-s12*s33-s23*s32+s13*s32+s12*s23-s13*s22) + m2*(-s21*s33+s11*s33+s23*s31-s13*s31-s11*s23+s13*s21) + m3*(s21*s32-s11*s32-s22*s31+s12*s31+s11*s22-s12*s21)
              
              x1 = x1 / det
              x2 = x2 / det
              x3 = x3 / det
              x4 = x4 / det
              
              if(min(x1,x2,x3,x4) >= 0.0 && max(x1,x2,x3,x4) <= 1.0)
              {
                sc <- sc + 1
              }
              else
              {
                snc <- snc + 1
              }
              
              w1 <- c(w1, x1)
              w2 <- c(w2, x2)
              w3 <- c(w3, x3)
              w4 <- c(w4, x4)
            }
          }
          
          df[nrow(df) + 1,] <- list(paste0(tracer[i],' ',tracer[j],' ',tracer[k]), w1[1], w2[1], w3[1], w4[1], (quantile(w1, 0.84)[[1]] - quantile(w1, 0.16)[[1]])/2, (quantile(w2, 0.84)[[1]] - quantile(w2, 0.16)[[1]])/2, (quantile(w3, 0.84)[[1]] - quantile(w3, 0.16)[[1]])/2, (quantile(w4, 0.84)[[1]] - quantile(w4, 0.16)[[1]])/2, sc/(sc+snc))
        }
      }
    }
    setTxtProgressBar(pb, i)}
  
  df <- transform(df, Dmax = pmax(Dw1, Dw2, Dw3, Dw4))
  df <- df[order(df$Dmax),]
  row.names(df) <- NULL
  cat('\n')
  cat('\n')
  return(df)
}
