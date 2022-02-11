#'CR values of each sediment mixture for three sources
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
cr_3s <- function(source, mixture, maxiter = 2000, seed = 123456)
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
		var <- sample(c(1:cols), 4)

		if(min(cr1[var]) < maxiter)
		{
			gof <- c()
			w1 <- c()
			w2 <- c()
			w3 <- c()
			rnd <- matrix(0, nrow = 3, ncol = 4)
			for (i in c(1:4))
			{
				j <- cols*2+1 # n column
				rnd[1, i] <- rt(1, source[1, j]) / sqrt(source[1, j])
				rnd[2, i] <- rt(1, source[2, j]) / sqrt(source[2, j])
				rnd[3, i] <- rt(1, source[3, j]) / sqrt(source[3, j])
			}
			for (i in c(1:4))
			{		
				j <- c(1:4)[-i]

				k <- j[1]
				l <- var[k]
				s11 <- source[1, l] + source[1, l+cols] * rnd[1, k]
				s12 <- source[2, l] + source[2, l+cols] * rnd[2, k]
				s13 <- source[3, l] + source[3, l+cols] * rnd[3, k]
				m1  <- mixture[l]

				k <- j[2]
				l <- var[k]
				s21 <- source[1, l] + source[1, l+cols] * rnd[1, k]
				s22 <- source[2, l] + source[2, l+cols] * rnd[2, k]
				s23 <- source[3, l] + source[3, l+cols] * rnd[3, k]
				m2  <- mixture[l]

				k <- j[3]
				l <- var[k]
				s31 <- source[1, l] + source[1, l+cols] * rnd[1, k]
				s32 <- source[2, l] + source[2, l+cols] * rnd[2, k]
				s33 <- source[3, l] + source[3, l+cols] * rnd[3, k]
				m3  <- mixture[l]

				m11 <- s33*s33-2.0*s31*s33+s31*s31+s23*s23-2.0*s21*s23+s21*s21+s13*s13-2.0*s11*s13+s11*s11
				m12 <- s33*s33-s32*s33-s31*s33+s31*s32+s23*s23-s22*s23-s21*s23+s21*s22+s13*s13-s12*s13-s11*s13+s11*s12
				m21 <- s33*s33-s32*s33-s31*s33+s31*s32+s23*s23-s22*s23-s21*s23+s21*s22+s13*s13-s12*s13-s11*s13+s11*s12
				m22 <- s33*s33-2.0*s32*s33+s32*s32+s23*s23-2.0*s22*s23+s22*s22+s13*s13-2.0*s12*s13+s12*s12

				v1 <- -(s33*m3-s31*m3+s23*m2-s21*m2+s13*m1-s11*m1-s33*s33+s31*s33-s23*s23+s21*s23-s13*s13+s11*s13)
				v2 <- -(s33*m3-s32*m3+s23*m2-s22*m2+s13*m1-s12*m1-s33*s33+s32*s33-s23*s23+s22*s23-s13*s13+s12*s13)

				det <- m11*m22-m12*m21

				if(det != 0.0)
				{
					n11 <- ( m22) / det
					n12 <- (-m12) / det
					n21 <- (-m21) / det
					n22 <- ( m11) / det

					x1 <- n11 * v1 + n12 * v2
					x2 <- n21 * v1 + n22 * v2
					x3 <- 1.0 - x1 - x2

					o1 <- s11 * x1 + s12 * x2 + s13 * x3
					o2 <- s21 * x1 + s22 * x2 + s23 * x3
					o3 <- s31 * x1 + s32 * x2 + s33 * x3

					err <- (o1-m1)*(o1-m1)+(o2-m2)*(o2-m2)+(o3-m3)*(o3-m3)
				}

				gof <- c(gof, err)
				w1 <- c(w1, x1)
				w2 <- c(w2, x2)
				w3 <- c(w3, x3)
			}

			sgof <- sort(gof)
			if(sgof[4]-sgof[3] > 10.0*(sgof[3]-sgof[1]))
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
				sol <- c(w1[win], w2[win], w3[win])
				if(max(sol) <= 1.0 && min(sol) >= 0.0)
				{
					for (i in c(1:4))
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