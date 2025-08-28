#' Quartets function 
#'
#' Extract all the possible selections formed by sets of n ??? 1 tracers and solve them by using standard methods
#'
#' @param source Data frame containing the sediment sources from a dataset
#' @param mixture Data frame containing one of the dataset mixtures
#' @param iter Number of iteration for each tracer
#' @param seed Seed for the random number generator
#' 
#' @return Data frame containing all the possible quartets combination from your dataset, their system of equation solution, the consistency and the discriminant capacity. 
#'
#' @export
#'
quartets <- function(source, mixture, iter = 1000, seed = 123456)
{
  set.seed(seed)
  
  source <- data.matrix(source[-1])
  mixture <- data.matrix(mixture[-1])
  
  cols <- (ncol(source)-1)/2
  tracer <- colnames(source)[1:cols]
  
  n <- cols*2+1 # n column
  
  df <- data.frame(id=character(), w1=double(), w2=double(), w3=double(), w4=double(), w5=double(), Dw1=double(), Dw2=double(), Dw3=double(), Dw4=double(), Dw5=double(), cons=double())
  
  # Introduce the progress bar
  pb <- txtProgressBar(min = 0, max = cols, style = 3, width = 50, char = "=")
  
  for (i in c(1:cols))
  {
    for (i2 in c(i+1:cols))
    {
      for (i3 in c(i2+1:cols))
      {
		    for (i4 in c(i3+1:cols))
		    {
		      if(i<i2 & i2<i3 & i3<i4 & i4<=cols)
		      {
		        w1 <- c()
		        w2 <- c()
		        w3 <- c()
		        w4 <- c()
		        w5 <- c()
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
		            s51 <- source[5, i]
		            m1  <- mixture[i]
		            
		            s12 <- source[1, i2]
		            s22 <- source[2, i2]
		            s32 <- source[3, i2]
		            s42 <- source[4, i2]
		            s52 <- source[5, i2]
		            m2  <- mixture[i2]
		            
		            s13 <- source[1, i3]
		            s23 <- source[2, i3]
		            s33 <- source[3, i3]
		            s43 <- source[4, i3]
		            s53 <- source[5, i3]
		            m3  <- mixture[i3]
		            
		            s14 <- source[1, i4]
		            s24 <- source[2, i4]
		            s34 <- source[3, i4]
		            s44 <- source[4, i4]
		            s54 <- source[5, i4]
		            m4  <- mixture[i4]
		          }
		          else
		          {
		            s11 <- source[1, i] + source[1, i+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
		            s21 <- source[2, i] + source[2, i+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
		            s31 <- source[3, i] + source[3, i+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
		            s41 <- source[4, i] + source[4, i+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
		            s51 <- source[5, i] + source[5, i+cols] * rt(1, source[5, n]) / sqrt(source[5, n])
		            m1  <- mixture[i]
		            
		            s12 <- source[1, i2] + source[1, i2+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
		            s22 <- source[2, i2] + source[2, i2+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
		            s32 <- source[3, i2] + source[3, i2+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
		            s42 <- source[4, i2] + source[4, i2+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
		            s52 <- source[5, i2] + source[5, i2+cols] * rt(1, source[5, n]) / sqrt(source[5, n])
		            m2  <- mixture[i2]
		            
		            s13 <- source[1, i3] + source[1, i3+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
		            s23 <- source[2, i3] + source[2, i3+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
		            s33 <- source[3, i3] + source[3, i3+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
		            s43 <- source[4, i3] + source[4, i3+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
		            s53 <- source[5, i3] + source[5, i3+cols] * rt(1, source[5, n]) / sqrt(source[5, n])
		            m3  <- mixture[i3]

		            s14 <- source[1, i4] + source[1, i4+cols] * rt(1, source[1, n]) / sqrt(source[1, n])
		            s24 <- source[2, i4] + source[2, i4+cols] * rt(1, source[2, n]) / sqrt(source[2, n])
		            s34 <- source[3, i4] + source[3, i4+cols] * rt(1, source[3, n]) / sqrt(source[3, n])
		            s44 <- source[4, i4] + source[4, i4+cols] * rt(1, source[4, n]) / sqrt(source[4, n])
		            s54 <- source[5, i4] + source[5, i4+cols] * rt(1, source[5, n]) / sqrt(source[5, n])
		            m4  <- mixture[i4]
		          }

							# maxima script
							#matrix(
							# [1,1,1,1,1],
							# [s11,s21,s31,s41,s51], 
							# [s12,s22,s32,s42,s52], 
							# [s13,s23,s33,s43,s53], 
							# [s14,s24,s34,s44,s54]
							#);
							#M:%;
							#determinant(M);
							#invert(M)*determinant(M);
							#ratsimp(%);

		          det <- s21*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-s11*(s32*(s43*s54-s44*s53)-s42*(s33*s54-s34*s53)+(s33*s44-s34*s43)*s52)-s31*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+s11*(s22*(s43*s54-s44*s53)-s42*(s23*s54-s24*s53)+(s23*s44-s24*s43)*s52)+s31*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)-s21*(s12*(s43*s54-s44*s53)-s42*(s13*s54-s14*s53)+(s13*s44-s14*s43)*s52)+s41*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-s11*(s22*(s33*s54-s34*s53)-s32*(s23*s54-s24*s53)+(s23*s34-s24*s33)*s52)-s41*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+s21*(s12*(s33*s54-s34*s53)-s32*(s13*s54-s14*s53)+(s13*s34-s14*s33)*s52)+s41*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-s31*(s12*(s23*s54-s24*s53)-s22*(s13*s54-s14*s53)+(s13*s24-s14*s23)*s52)-(s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)*s51+(s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)*s51-(s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)*s51+(s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s51+s11*(s22*(s33*s44-s34*s43)-s32*(s23*s44-s24*s43)+(s23*s34-s24*s33)*s42)-s21*(s12*(s33*s44-s34*s43)-s32*(s13*s44-s14*s43)+(s13*s34-s14*s33)*s42)+s31*(s12*(s23*s44-s24*s43)-s22*(s13*s44-s14*s43)+(s13*s24-s14*s23)*s42)-(s12*(s23*s34-s24*s33)-s22*(s13*s34-s14*s33)+(s13*s24-s14*s23)*s32)*s41
	
		          if(det != 0.0)
		          {
		            x1 = ((s21*s32-s22*s31)*s43+(s23*s31-s21*s33)*s42+(s22*s33-s23*s32)*s41)*s54+((s22*s31-s21*s32)*s44+(s21*s34-s24*s31)*s42+(s24*s32-s22*s34)*s41)*s53+((s21*s33-s23*s31)*s44+(s24*s31-s21*s34)*s43+(s23*s34-s24*s33)*s41)*s52+((s23*s32-s22*s33)*s44+(s22*s34-s24*s32)*s43+(s24*s33-s23*s34)*s42)*s51 + m1*(((s22-s32)*s43+(s33-s23)*s42-s22*s33+s23*s32)*s54+((s32-s22)*s44+(s24-s34)*s42+s22*s34-s24*s32)*s53+((s23-s33)*s44+(s34-s24)*s43-s23*s34+s24*s33)*s52+(s22*s33-s23*s32)*s44+(s24*s32-s22*s34)*s43+(s23*s34-s24*s33)*s42) + m2*(((s31-s21)*s43+(s23-s33)*s41+s21*s33-s23*s31)*s54+((s21-s31)*s44+(s34-s24)*s41-s21*s34+s24*s31)*s53+((s33-s23)*s44+(s24-s34)*s43+s23*s34-s24*s33)*s51+(s23*s31-s21*s33)*s44+(s21*s34-s24*s31)*s43+(s24*s33-s23*s34)*s41) + m3*(((s21-s31)*s42+(s32-s22)*s41-s21*s32+s22*s31)*s54+((s31-s21)*s44+(s24-s34)*s41+s21*s34-s24*s31)*s52+((s22-s32)*s44+(s34-s24)*s42-s22*s34+s24*s32)*s51+(s21*s32-s22*s31)*s44+(s24*s31-s21*s34)*s42+(s22*s34-s24*s32)*s41) + m4*(((s31-s21)*s42+(s22-s32)*s41+s21*s32-s22*s31)*s53+((s21-s31)*s43+(s33-s23)*s41-s21*s33+s23*s31)*s52+((s32-s22)*s43+(s23-s33)*s42+s22*s33-s23*s32)*s51+(s22*s31-s21*s32)*s43+(s21*s33-s23*s31)*s42+(s23*s32-s22*s33)*s41)
		            
		            x2 = ((s12*s31-s11*s32)*s43+(s11*s33-s13*s31)*s42+(s13*s32-s12*s33)*s41)*s54+((s11*s32-s12*s31)*s44+(s14*s31-s11*s34)*s42+(s12*s34-s14*s32)*s41)*s53+((s13*s31-s11*s33)*s44+(s11*s34-s14*s31)*s43+(s14*s33-s13*s34)*s41)*s52+((s12*s33-s13*s32)*s44+(s14*s32-s12*s34)*s43+(s13*s34-s14*s33)*s42)*s51 + m1*(((s32-s12)*s43+(s13-s33)*s42+s12*s33-s13*s32)*s54+((s12-s32)*s44+(s34-s14)*s42-s12*s34+s14*s32)*s53+((s33-s13)*s44+(s14-s34)*s43+s13*s34-s14*s33)*s52+(s13*s32-s12*s33)*s44+(s12*s34-s14*s32)*s43+(s14*s33-s13*s34)*s42) + m2*(((s11-s31)*s43+(s33-s13)*s41-s11*s33+s13*s31)*s54+((s31-s11)*s44+(s14-s34)*s41+s11*s34-s14*s31)*s53+((s13-s33)*s44+(s34-s14)*s43-s13*s34+s14*s33)*s51+(s11*s33-s13*s31)*s44+(s14*s31-s11*s34)*s43+(s13*s34-s14*s33)*s41) + m3*(((s31-s11)*s42+(s12-s32)*s41+s11*s32-s12*s31)*s54+((s11-s31)*s44+(s34-s14)*s41-s11*s34+s14*s31)*s52+((s32-s12)*s44+(s14-s34)*s42+s12*s34-s14*s32)*s51+(s12*s31-s11*s32)*s44+(s11*s34-s14*s31)*s42+(s14*s32-s12*s34)*s41) + m4*(((s11-s31)*s42+(s32-s12)*s41-s11*s32+s12*s31)*s53+((s31-s11)*s43+(s13-s33)*s41+s11*s33-s13*s31)*s52+((s12-s32)*s43+(s33-s13)*s42-s12*s33+s13*s32)*s51+(s11*s32-s12*s31)*s43+(s13*s31-s11*s33)*s42+(s12*s33-s13*s32)*s41)
		            
		            x3 = ((s11*s22-s12*s21)*s43+(s13*s21-s11*s23)*s42+(s12*s23-s13*s22)*s41)*s54+((s12*s21-s11*s22)*s44+(s11*s24-s14*s21)*s42+(s14*s22-s12*s24)*s41)*s53+((s11*s23-s13*s21)*s44+(s14*s21-s11*s24)*s43+(s13*s24-s14*s23)*s41)*s52+((s13*s22-s12*s23)*s44+(s12*s24-s14*s22)*s43+(s14*s23-s13*s24)*s42)*s51 + m1*(((s12-s22)*s43+(s23-s13)*s42-s12*s23+s13*s22)*s54+((s22-s12)*s44+(s14-s24)*s42+s12*s24-s14*s22)*s53+((s13-s23)*s44+(s24-s14)*s43-s13*s24+s14*s23)*s52+(s12*s23-s13*s22)*s44+(s14*s22-s12*s24)*s43+(s13*s24-s14*s23)*s42) + m2*(((s21-s11)*s43+(s13-s23)*s41+s11*s23-s13*s21)*s54+((s11-s21)*s44+(s24-s14)*s41-s11*s24+s14*s21)*s53+((s23-s13)*s44+(s14-s24)*s43+s13*s24-s14*s23)*s51+(s13*s21-s11*s23)*s44+(s11*s24-s14*s21)*s43+(s14*s23-s13*s24)*s41) + m3*(((s11-s21)*s42+(s22-s12)*s41-s11*s22+s12*s21)*s54+((s21-s11)*s44+(s14-s24)*s41+s11*s24-s14*s21)*s52+((s12-s22)*s44+(s24-s14)*s42-s12*s24+s14*s22)*s51+(s11*s22-s12*s21)*s44+(s14*s21-s11*s24)*s42+(s12*s24-s14*s22)*s41) + m4*(((s21-s11)*s42+(s12-s22)*s41+s11*s22-s12*s21)*s53+((s11-s21)*s43+(s23-s13)*s41-s11*s23+s13*s21)*s52+((s22-s12)*s43+(s13-s23)*s42+s12*s23-s13*s22)*s51+(s12*s21-s11*s22)*s43+(s11*s23-s13*s21)*s42+(s13*s22-s12*s23)*s41)
		            
		            x4 = ((s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)*s54+((s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31)*s53+((s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31)*s52+((s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32)*s51 + m1*(((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s54+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s53+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s52+(s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32) + m2*(((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s54+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s53+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s51+(s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31) + m3*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s54+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s52+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s51+(s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31) + m4*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s53+((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s52+((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s51+(s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)
		            
		            x5 = ((s11*s22-s12*s21)*s33+(s13*s21-s11*s23)*s32+(s12*s23-s13*s22)*s31)*s44+((s12*s21-s11*s22)*s34+(s11*s24-s14*s21)*s32+(s14*s22-s12*s24)*s31)*s43+((s11*s23-s13*s21)*s34+(s14*s21-s11*s24)*s33+(s13*s24-s14*s23)*s31)*s42+((s13*s22-s12*s23)*s34+(s12*s24-s14*s22)*s33+(s14*s23-s13*s24)*s32)*s41 + m1*(((s12-s22)*s33+(s23-s13)*s32-s12*s23+s13*s22)*s44+((s22-s12)*s34+(s14-s24)*s32+s12*s24-s14*s22)*s43+((s13-s23)*s34+(s24-s14)*s33-s13*s24+s14*s23)*s42+(s12*s23-s13*s22)*s34+(s14*s22-s12*s24)*s33+(s13*s24-s14*s23)*s32) + m2*(((s21-s11)*s33+(s13-s23)*s31+s11*s23-s13*s21)*s44+((s11-s21)*s34+(s24-s14)*s31-s11*s24+s14*s21)*s43+((s23-s13)*s34+(s14-s24)*s33+s13*s24-s14*s23)*s41+(s13*s21-s11*s23)*s34+(s11*s24-s14*s21)*s33+(s14*s23-s13*s24)*s31) + m3*(((s11-s21)*s32+(s22-s12)*s31-s11*s22+s12*s21)*s44+((s21-s11)*s34+(s14-s24)*s31+s11*s24-s14*s21)*s42+((s12-s22)*s34+(s24-s14)*s32-s12*s24+s14*s22)*s41+(s11*s22-s12*s21)*s34+(s14*s21-s11*s24)*s32+(s12*s24-s14*s22)*s31) + m4*(((s21-s11)*s32+(s12-s22)*s31+s11*s22-s12*s21)*s43+((s11-s21)*s33+(s23-s13)*s31-s11*s23+s13*s21)*s42+((s22-s12)*s33+(s13-s23)*s32+s12*s23-s13*s22)*s41+(s12*s21-s11*s22)*s33+(s11*s23-s13*s21)*s32+(s13*s22-s12*s23)*s31)
		            
		            x1 = x1 / det
		            x2 = x2 / det
		            x3 = x3 / det
		            x4 = x4 / det
		            x5 = x5 / det
		            
		            if(min(x1,x2,x3,x4,x5) >= 0.0 && max(x1,x2,x3,x4,x5) <= 1.0)
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
		            w5 <- c(w5, x5)
		          }
		        }
		        
		       	# df[nrow(df) + 1,] <- list(paste0(tracer[i],' ',tracer[i1],' ',tracer[i2],' ',tracer[i3]), w1[1], w2[1], w3[1], w4[1], w5[1], (quantile(w1, 0.84)[[1]] - quantile(w1, 0.16)[[1]])/2, (quantile(w2, 0.84)[[1]] - quantile(w2, 0.16)[[1]])/2, (quantile(w3, 0.84)[[1]] - quantile(w3, 0.16)[[1]])/2, (quantile(w4, 0.84)[[1]] - quantile(w4, 0.16)[[1]])/2,  (quantile(w5, 0.84)[[1]] - quantile(w5, 0.16)[[1]])/2, sc/(sc+snc))
		       	df[nrow(df) + 1,] <- list(paste0(tracer[i],' ',tracer[i2],' ',tracer[i3],' ',tracer[i4]), w1[1], w2[1], w3[1], w4[1], w5[1], (quantile(w1, 0.84)[[1]] - quantile(w1, 0.16)[[1]])/2, (quantile(w2, 0.84)[[1]] - quantile(w2, 0.16)[[1]])/2, (quantile(w3, 0.84)[[1]] - quantile(w3, 0.16)[[1]])/2, (quantile(w4, 0.84)[[1]] - quantile(w4, 0.16)[[1]])/2,  (quantile(w5, 0.84)[[1]] - quantile(w5, 0.16)[[1]])/2, sc/(sc+snc))
		      }
        }
      }
    }
    setTxtProgressBar(pb, i)}
  
  df <- transform(df, Dmax = pmax(Dw1, Dw2, Dw3, Dw4, Dw5))
  df <- df[order(df$Dmax),]
  row.names(df) <- NULL
  cat('\n')
  cat('\n')
  return(df)
}

