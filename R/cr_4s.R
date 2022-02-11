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
cr_4s <- function(source, mixture, maxiter = 2000, seed = 123456)
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
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = maxiter, style = 3, width = 50, char = "=")
  
  while (min(cr1) < maxiter)
  {
    var <- sample(c(1:cols), 5)
    
    if(min(cr1[var]) < maxiter)
    {
      gof <- c()
      w1 <- c()
      w2 <- c()
      w3 <- c()
      w4 <- c()
      rnd <- matrix(0, nrow = 4, ncol = 5)
      for (i in c(1:5))
      {
        j <- cols*2+1 # n column
        rnd[1, i] <- rt(1, source[1, j]) / sqrt(source[1, j])
        rnd[2, i] <- rt(1, source[2, j]) / sqrt(source[2, j])
        rnd[3, i] <- rt(1, source[3, j]) / sqrt(source[3, j])
        rnd[4, i] <- rt(1, source[4, j]) / sqrt(source[4, j])
      }
      for (i in c(1:5))
      {		
        j <- c(1:5)[-i]
        
        k <- j[1]
        l <- var[k]
        s11 <- source[1, l] + source[1, l+cols] * rnd[1, k]
        s12 <- source[2, l] + source[2, l+cols] * rnd[2, k]
        s13 <- source[3, l] + source[3, l+cols] * rnd[3, k]
        s14 <- source[4, l] + source[4, l+cols] * rnd[4, k]
        m1  <- mixture[l]
        
        k <- j[2]
        l <- var[k]
        s21 <- source[1, l] + source[1, l+cols] * rnd[1, k]
        s22 <- source[2, l] + source[2, l+cols] * rnd[2, k]
        s23 <- source[3, l] + source[3, l+cols] * rnd[3, k]
        s24 <- source[4, l] + source[4, l+cols] * rnd[4, k]
        m2  <- mixture[l]
        
        k <- j[3]
        l <- var[k]
        s31 <- source[1, l] + source[1, l+cols] * rnd[1, k]
        s32 <- source[2, l] + source[2, l+cols] * rnd[2, k]
        s33 <- source[3, l] + source[3, l+cols] * rnd[3, k]
        s34 <- source[4, l] + source[4, l+cols] * rnd[4, k]
        m3  <- mixture[l]
        
        k <- j[4]
        l <- var[k]
        s41 <- source[1, l] + source[1, l+cols] * rnd[1, k]
        s42 <- source[2, l] + source[2, l+cols] * rnd[2, k]
        s43 <- source[3, l] + source[3, l+cols] * rnd[3, k]
        s44 <- source[4, l] + source[4, l+cols] * rnd[4, k]
        m4  <- mixture[l]
        
        m11 <- s44*s44-2.0*s41*s44+s41*s41+s34*s34-2.0*s31*s34+s31*s31+s24*s24-2.0*s21*s24+s21*s21+s14*s14-2.0*s11*s14+s11*s11
        m12 <- s44*s44-s42*s44-s41*s44+s41*s42+s34*s34-s32*s34-s31*s34+s31*s32+s24*s24-s22*s24-s21*s24+s21*s22+s14*s14-s12*s14-s11*s14+s11*s12
        m13 <- s44*s44-s43*s44-s41*s44+s41*s43+s34*s34-s33*s34-s31*s34+s31*s33+s24*s24-s23*s24-s21*s24+s21*s23+s14*s14-s13*s14-s11*s14+s11*s13
        m21 <- s44*s44-s42*s44-s41*s44+s41*s42+s34*s34-s32*s34-s31*s34+s31*s32+s24*s24-s22*s24-s21*s24+s21*s22+s14*s14-s12*s14-s11*s14+s11*s12
        m22 <- s44*s44-2.0*s42*s44+s42*s42+s34*s34-2.0*s32*s34+s32*s32+s24*s24-2.0*s22*s24+s22*s22+s14*s14-2.0*s12*s14+s12*s12
        m23 <- s44*s44-s43*s44-s42*s44+s42*s43+s34*s34-s33*s34-s32*s34+s32*s33+s24*s24-s23*s24-s22*s24+s22*s23+s14*s14-s13*s14-s12*s14+s12*s13
        m31 <- s44*s44-s43*s44-s41*s44+s41*s43+s34*s34-s33*s34-s31*s34+s31*s33+s24*s24-s23*s24-s21*s24+s21*s23+s14*s14-s13*s14-s11*s14+s11*s13
        m32 <- s44*s44-s43*s44-s42*s44+s42*s43+s34*s34-s33*s34-s32*s34+s32*s33+s24*s24-s23*s24-s22*s24+s22*s23+s14*s14-s13*s14-s12*s14+s12*s13
        m33 <- s44*s44-2.0*s43*s44+s43*s43+s34*s34-2.0*s33*s34+s33*s33+s24*s24-2.0*s23*s24+s23*s23+s14*s14-2.0*s13*s14+s13*s13
        
        v1 <- -(s44*m4-s41*m4+s34*m3-s31*m3+s24*m2-s21*m2+s14*m1-s11*m1-s44*s44+s41*s44-s34*s34+s31*s34-s24*s24+s21*s24-s14*s14+s11*s14)
        v2 <- -(s44*m4-s42*m4+s34*m3-s32*m3+s24*m2-s22*m2+s14*m1-s12*m1-s44*s44+s42*s44-s34*s34+s32*s34-s24*s24+s22*s24-s14*s14+s12*s14)
        v3 <- -(s44*m4-s43*m4+s34*m3-s33*m3+s24*m2-s23*m2+s14*m1-s13*m1-s44*s44+s43*s44-s34*s34+s33*s34-s24*s24+s23*s24-s14*s14+s13*s14)
        
        det <- (m11*m22-m12*m21)*m33+(m13*m21-m11*m23)*m32+(m12*m23-m13*m22)*m31
        
        if(det != 0.0)
        {
          n11 <- (m22*m33-m23*m32)/det
          n12 <- (m13*m32-m12*m33)/det
          n13 <- (m12*m23-m13*m22)/det
          n21 <- (m23*m31-m21*m33)/det
          n22 <- (m11*m33-m13*m31)/det
          n23 <- (m13*m21-m11*m23)/det
          n31 <- (m21*m32-m22*m31)/det
          n32 <- (m12*m31-m11*m32)/det
          n33 <- (m11*m22-m12*m21)/det
          
          x1 <- n11 * v1 + n12 * v2 + n13 * v3
          x2 <- n21 * v1 + n22 * v2 + n23 * v3
          x3 <- n31 * v1 + n32 * v2 + n33 * v3
          x4 <- 1.0 - x1 - x2 - x3
          
          o1 <- s11 * x1 + s12 * x2 + s13 * x3 + s14 * x4
          o2 <- s21 * x1 + s22 * x2 + s23 * x3 + s24 * x4
          o3 <- s31 * x1 + s32 * x2 + s33 * x3 + s34 * x4
          o4 <- s41 * x1 + s42 * x2 + s43 * x3 + s44 * x4
          
          err <- (o1-m1)*(o1-m1)+(o2-m2)*(o2-m2)+(o3-m3)*(o3-m3)+(o4-m4)*(o4-m4)
        }
        
        gof <- c(gof, err)
        w1 <- c(w1, x1)
        w2 <- c(w2, x2)
        w3 <- c(w3, x3)
        w4 <- c(w4, x4)
      }
      
      sgof <- sort(gof)
      if(sgof[5]-sgof[4] > 10.0*(sgof[4]-sgof[1]))
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
        sol <- c(w1[win], w2[win], w3[win], w4[win])
        if(max(sol) <= 1.0 && min(sol) >= 0.0)
        {
          for (i in c(1:5))
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