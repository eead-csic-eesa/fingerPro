#' CR values of each sediment mixture
#'
#' Compute the consensus ranking method for the selected dataset.
#'
#' @param data Data frame containing sediment source and mixtures
#' @param maxiter Number of iteraction for each tracer
#'
#' @return Data frame containing the CR score for each tracer.
#'
#' @export
#'
CR_Method <- function(data, niter = 100L, seed = 123456L){

# normalise 
cols <- colnames(data)[-1][-1]
for (col in cols)
{
	mx <- max(data[col])
	mn <- min(data[col])
	data[col] <- ( data[col] - mn ) / ( mx - mn )
}

cr1 <- c()
cr2 <- c()
for (v in colnames(data))
{
	cr1[v] = 0
	cr2[v] = 0
}

allsources <- inputSource(data)
allmixtures <- inputSample(data)

maxiter <- niter # at least 1000, however it can be tried with 100 to check the general trend
for (iter in c(1:100000))
{
	var <- sample(colnames(data)[-1][-1], 4)

	run <- FALSE
	for (i in c(1:4))
	{
		if(cr1[var[i]] < maxiter)
		{
			run <- TRUE
		}
	}

	if(run)
	{
		gof <- c()
		w1 <- c()
		w2 <- c()
		w3 <- c()
		for (i in c(1:4))
		{
			sources  <-  allsources[c("id", var[-i], paste0("D", var[-i]), "n")]
			mixtures <- allmixtures[c("id", var[-i])]

			s11 <- sources[[2]][1] + sources[[5]][1] * rt(1, sources[[8]][1] - 1) / sqrt(sources[[8]][1])
			s12 <- sources[[2]][2] + sources[[5]][2] * rt(1, sources[[8]][2] - 1) / sqrt(sources[[8]][2])
			s13 <- sources[[2]][3] + sources[[5]][3] * rt(1, sources[[8]][3] - 1) / sqrt(sources[[8]][3])
			m1  <- mixtures[[2]]

			s21 <- sources[[3]][1] + sources[[6]][1] * rt(1, sources[[8]][1] - 1) / sqrt(sources[[8]][1])
			s22 <- sources[[3]][2] + sources[[6]][2] * rt(1, sources[[8]][2] - 1) / sqrt(sources[[8]][2])
			s23 <- sources[[3]][3] + sources[[6]][3] * rt(1, sources[[8]][3] - 1) / sqrt(sources[[8]][3])
			m2  <- mixtures[[3]]

			s31 <- sources[[4]][1] + sources[[7]][1] * rt(1, sources[[8]][1] - 1) / sqrt(sources[[8]][1])
			s32 <- sources[[4]][2] + sources[[7]][2] * rt(1, sources[[8]][2] - 1) / sqrt(sources[[8]][2])
			s33 <- sources[[4]][3] + sources[[7]][3] * rt(1, sources[[8]][3] - 1) / sqrt(sources[[8]][3])
			m3  <- mixtures[[4]]

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
		win <- which.min(gof)

		sol <- c(w1[win], w2[win], w3[win])
		if(max(sol) <= 1.0 && min(sol) >= 0.0)
		{
			for (i in c(1:4))
			{
				if(cr1[var[i]] < maxiter)
				{
					cr1[var[i]] = cr1[var[i]] + 1

					if(i == win)
					{
						cr2[var[i]] = cr2[var[i]] + 1
					}
				}
			}
		}
	}
}
CR_val <- 100-((cr2*100)/cr1)

CR_values <- cbind.data.frame(cr1,cr2,CR_val)
colnames(CR_values) <- c("Total debates", "Lost debates", "CR_Values")
CR_values <- as.data.frame(t(CR_values[3:nrow(CR_values),]))

print(CR_values)
return(CR_values)
}

cr_4s <- function(source, mixture, maxiter=10, seed=123456)
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
	}

	score <- 100-((cr2*100)/cr1)
	cr <- data.frame(tracer, score)
	cr <- cr[rev(order(cr$score)),]
	row.names(cr) <- NULL
	
	return(cr)
}


cr_3s <- function(source, mixture, maxiter=10, seed=123456)
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
	}

	score <- 100-((cr2*100)/cr1)
	cr <- data.frame(tracer, score)
	cr <- cr[rev(order(cr$score)),]
	row.names(cr) <- NULL
	
	return(cr)
}

