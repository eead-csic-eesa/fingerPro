cts_4s <- function(source, mixture, sol)
{
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
	
	err <- c()
	for (col in c(1:cols))
	{
		err[col] <- abs(source[1,col] * sol[1] + source[2,col] * sol[2] + source[3,col] * sol[3] + source[4,col] * sol[4] -  mixture[1,col])
	}
	
	return(data.frame(tracer, err))
}

cts_3s <- function(source, mixture, sol)
{
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
	
	err <- c()
	for (col in c(1:cols))
	{
		err[col] <- abs(source[1,col] * sol[1] + source[2,col] * sol[2] + source[3,col] * sol[3] -  mixture[1,col])
	}
	
	return(data.frame(tracer, err))
}