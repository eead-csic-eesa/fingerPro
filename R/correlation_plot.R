#' Correlation matrix chart
#'
#' The function displays a correlation matrix of each of the properties divided by the different sources to help the user in the decision.
#'
#' @param data Data frame containing source and mixtures data
#' @param columns Numeric vector containing the index of the columns in the chart (the first column refers to the grouping variable)
#' @param mixtures Boolean to include or exclude the mixture samples in the chart
#' @param nmixtures Number of mixtures in the dataset
#' @param colors Vector of colors to use for the scatterplot
#'
#' @export
#'
correlationPlot <- function(data, columns = c(1:ncol(data)-1), mixtures = FALSE, nmixtures = 1, colors = NULL)  {
  # reorder the groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  data1 <- data
  colnames(data1)[2] <- "Sources"
  if (mixtures) {
    # Here three mixture points are created to highlight the mixture samples in the graph
    if (nmixtures < 3) {
      data1[nrow(data1) + 1,] <- data1[nrow(data1),]
      data1[nrow(data1) + 1,] <- data1[nrow(data1),]
    }
  } else {
    data1 <- data1[which(!data1[,2] == data1[nrow(data1),2]), ]
  }
  
  if (!is.null(colors)) {
    p <- ggpairs(data = data1, columns = columns + 1, upper = list(continuous = "cor"), 
                 lower = list(combo = "facetdensity"), title = "Correlation Matrix", 
                 mapping = aes(color = Sources, alpha = 0.25), legend = c(1, 1)) +
      theme_bw() +
      theme(legend.position = "bottom") +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors)
  } else {
    p <- ggpairs(data = data1, columns = columns + 1, upper = list(continuous = "cor"), 
                 lower = list(combo = "facetdensity"), title = "Correlation Matrix", 
                 mapping = aes(color = Sources, alpha = 0.25), legend = c(1, 1)) +
      theme_bw() +
      theme(legend.position = "bottom")
  }
  
  return(p)
}
