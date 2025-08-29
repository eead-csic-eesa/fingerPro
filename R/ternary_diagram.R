#' @title Visualize individual tracer analysis as ternary diagrams
#'
#' @description This function creates ternary diagrams to visualize the results of the individual tracer analysis. Each ternary diagram represents the predicted apportionments for a specific tracer.
#'
#' @param data A data frame containing the results from the individual tracer analysis function.
#' @param tracers A vector specifying the indices of the tracers to be displayed.
#' @param rows An integer specifying the number of rows in the grid of ternary diagrams.
#' @param cols An integer specifying the number of columns in the grid of ternary diagrams.
#' @param solution A vector containing an optional reference solution for visual comparison.
#'
#' @return A grid of ternary diagrams, each representing the predicted apportionments for a specific tracer. If there are three sources, the function generates one ternary triangle for each tracer. If there are four sources, the function generates six triangles for each tracer. The six triangles represent the following source combinations at their vertices:
#' 1. (S1, S2, S3+S4)
#' 2. (S2, S3, S1+S4)
#' 3. (S3, S4, S1+S2)
#' 4. (S4, S1, S2+S3)
#' 5. (S1, S3, S2+S4)
#' 6. (S2, S4, S1+S3)
#'
#' @export
ternary_diagram <- function(data, tracers = c(1:2), rows = 1, cols = 2, solution = NA) {
  
  source_n <- nrow(inputSource(data[[length(data) - 1]]))
  
  if (source_n == 3) {
    t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
    t_names <- gsub("^mean_", "", t_names)
    
    plots <- list()
    par(mar = c(0, 0.1, 3, 0.1), mfcol = c(rows, cols))
    for (i2 in tracers) {
      result <- data[[i2]]
      labels <- colnames(result)
      labels <- gsub("^w.", "", labels)
      colnames(result) <- c("ID", "GOF", "w.S1", "w.S2", "w.S3")
      
      x <- result$w.S1
      y <- result$w.S2
      z <- result$w.S3
      
      test <- as.data.frame(cbind(x, y, z))

			if (length(solution)==3) {
				x_sol <- solution[1]
				y_sol <- solution[2]
				z_sol <- solution[3]
				test_sol <- as.data.frame(cbind(x_sol, y_sol, z_sol))
			}

      plots[[i2]] <- Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1) #, alab = labels[3], blab = labels[4], clab = labels[5])
      Ternary::TernaryPoints(test, col = scales::alpha('blue', 0.5), cex = 0.1)
			if (length(solution)==3) {
				Ternary::TernaryPoints(test_sol, col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
 			}
      title(t_names[i2], cex.main = 3.5)
		}
	} else if (source_n == 4) {
    t_names <- colnames(data[[length(data) - 1]][3:ncol(data[[length(data) - 1]])])
    t_names <- gsub("^mean_", "", t_names)
        
    plots <- list()
    par(mar = c(0, 0.1, 3, 0.1), mfcol = c(6, cols))
    for (i2 in tracers) {
      result <- data[[i2]]
      labels <- colnames(result)
      labels <- gsub("^w.", "", labels)
      colnames(result) <- c("ID", "GOF", "w.S1", "w.S2", "w.S3", "w.S4")
      
      # Original triangles
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
      
      test <- as.data.frame(cbind(x, y, z))
      test1 <- as.data.frame(cbind(x1, y1, z1))
      test2 <- as.data.frame(cbind(x2, y2, z2))
      test3 <- as.data.frame(cbind(x3, y3, z3))
      test4 <- as.data.frame(cbind(x4, y4, z4))
      test5 <- as.data.frame(cbind(x5, y5, z5))
      
      if (!is.na(solution)) {
		    # Virtual sample
		    x_sol <- solution[1]
		    y_sol <- solution[2]
		    z_sol <- solution[3] + solution[4]
		    test_sol <- as.data.frame(cbind(x_sol, y_sol, z_sol))
		    x_sol1 <- solution[2]
		    y_sol1 <- solution[3]
		    z_sol1 <- solution[1] + solution[4]
		    test_sol1 <- as.data.frame(cbind(x_sol1, y_sol1, z_sol1))
		    x_sol2 <- solution[3]
		    y_sol2 <- solution[4]
		    z_sol2 <- solution[1] + solution[2]
		    test_sol2 <- as.data.frame(cbind(x_sol2, y_sol2, z_sol2))
		    x_sol3 <- solution[4]
		    y_sol3 <- solution[1]
		    z_sol3 <- solution[2] + solution[3]
		    test_sol3 <- as.data.frame(cbind(x_sol3, y_sol3, z_sol3))
		    x_sol4 <- solution[1]
		    y_sol4 <- solution[3]
		    z_sol4 <- solution[2] + solution[4]
		    test_sol4 <- as.data.frame(cbind(x_sol4, y_sol4, z_sol4))
		    x_sol5 <- solution[2]
		    y_sol5 <- solution[4]
		    z_sol5 <- solution[1] + solution[3]
		    test_sol5 <- as.data.frame(cbind(x_sol5, y_sol5, z_sol5))
      }
      
      # Plotting
      # par(mfcol=c(6, cols))  # Set up multi-plot layout
      # Plot original triangles
      plots[[i2]] <- Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      title(t_names[i2], cex.main = 3.5)
      
      Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test1, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol1,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test2, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol2,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test3, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol3,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
       Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test4, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol4,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
      Ternary::TernaryPlot(xlim = c(-0.6, 0.6), tip.cex = 1, grid.lines = 4, grid.minor.lines = 1)
      Ternary::TernaryPoints(test5, col = scales::alpha('blue', 0.5), cex = 0.1)
      if (!is.na(solution)) {
        Ternary::TernaryPoints(test_sol5,  col = scales::alpha('red', 0.9), cex = 0.8, pch = 23, bg = "black", lwd = 2)
      }
    }
  } else {
  	stop(paste0("Error: ternary_diagram is not implemented for ", source_n," sources."))
  }
}
