#' Ternary diagrams with the CI function results
#'
#' Display the CI prediction for the selected number 
#'
#' @param data Data frame containing the results from the CI_Method function
#' @param tracers Number of tracers to display
#' @param Means Number of iteration for each tracer
#' @param n_col Number of columns in which display the tracers
#' @param n_row Number of rows in which display the tracers
#'
#' @return Ternary diagrams from the individual tracer predictions
#'
#' @export
#'
Ternary_diagram <- function(data, tracers= 1:2, n_row = 1, n_col = 2, Means = T){
  
  sources <- nrow(inputSource(data[[length(data)-1]]))
  
  if (sources == 3) {
    if (Means == T) {
      t_names <- colnames(data[[length(data)-1]][3:(ncol(data[[length(data)-1]]) - 3)])
    } else {
      t_names <- colnames(data[[length(data)-1]][3:ncol(data[[length(data)-1]])])
    }
    
  plots <- list()
  par(mar=c(0, 0.1, 3, 0.1),mfcol=c(n_row, n_col))
  for (i2 in tracers) {
    result <- data[[i2]]
    
    x <- result$w.S1
    y <- result$w.S2
    z <- result$w.S3
    
    test <- as.data.frame(cbind(x,y,z))
    plots [[i2]] <- TernaryPlot(atip = "S1", btip = "S2", ctip = "S3",xlim = c(-0.75,0.75)) 
      TernaryPoints(test, col = alpha('blue', 0.5), cex = 0.1) +
      title(t_names[i2], cex.main = 3.5)
    }
  } else if (sources == 4){
    
    if (Means == T) {
      t_names <- colnames(data[[length(data)-1]][3:(ncol(data[[length(data)-1]]) - 3)])
    } else {
      t_names <- colnames(data[[length(data)-1]][3:ncol(data[[length(data)-1]])])
    }
    par(mar=c(0, 0.1, 3, 0.1),mfcol=c(6,n_col))
    for (i2 in tracers) {
      result <- data[[i2]]
      
      x <- result$w.S1
      y <- result$w.S2
      z <- result$w.S3 + result$w.S4
      
      x1 <- result$w.S2
      y1 <- result$w.S3
      z1 <- result$w.S1 + result$w.S4
      
      x2 <- result$w.S3
      y2 <- result$w.S4
      z2 <- result$w.S1 + result$w.S2
      
      x3 <- result$w.S4
      y3 <- result$w.S1
      z3 <- result$w.S2 + result$w.S3
      
      x4 <- result$w.S1
      y4 <- result$w.S3
      z4 <- result$w.S2 + result$w.S4 
      
      x5 <- result$w.S2
      y5 <- result$w.S4
      z5 <- result$w.S1 + result$w.S3
      
      
      test <- as.data.frame(cbind(x,y,z))
      test1 <- as.data.frame(cbind(x1,y1,z1))
      test2 <- as.data.frame(cbind(x2,y2,z2))
      test3 <- as.data.frame(cbind(x3,y3,z3))
      test4 <- as.data.frame(cbind(x4,y4,z4))
      test5 <- as.data.frame(cbind(x5,y5,z5))

      
      TernaryPlot(atip="S1", btip = "S2", ctip = "S3 + S4",xlim = c(-0.75,0.75),tip.cex = 1) 
        TernaryPoints(test2, col=alpha('blue',0.5), cex=0.1) + 
        title(t_names[i2], cex.main=3.5, )
      TernaryPlot(atip="S2", btip = "S3", ctip = "S1 + S4",xlim = c(-0.75,0.75)) 
        TernaryPoints(test1, col=alpha('blue',0.5), cex=0.1) 
      TernaryPlot(atip="S3", btip = "S4", ctip = "S1 + S2",xlim = c(-0.75,0.75)) 
        TernaryPoints(test2, col=alpha('blue',0.5), cex=0.1) 
      TernaryPlot(atip="S4", btip = "S1", ctip = "S3 + S2",xlim = c(-0.75,0.75)) 
        TernaryPoints(test3, col=alpha('blue',0.5), cex=0.1) 
      TernaryPlot(atip="S1", btip = "S3", ctip = "S4 + S2",xlim = c(-0.75,0.75)) 
        TernaryPoints(test4, col=alpha('blue',0.5), cex=0.1) 
      TernaryPlot(atip="S2", btip = "S4", ctip = "S1 + S3",xlim = c(-0.75,0.75)) 
        TernaryPoints(test5, col=alpha('blue',0.5), cex=0.1) 
      
    }
  }
}
