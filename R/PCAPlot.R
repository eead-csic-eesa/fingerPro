#' Principal component analysis chart
#'
#' The function performs a principal components analysis on the given data matrix and displays a biplot using vqv.ggbiplot package of the results for each different source to help the user in the decision.
#'
#' @param data Data frame containing source and mixtures data
#' @param components Numeric vector containing the index of the two principal components in the chart
#'
#' @export
#' 

if(getRversion() >= "2.15.1")  utils::globalVariables(c("groups"))

PCAPlot <- function(data, components = c(1:2)) {

   # reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  data_PCA <- data[c(3:ncol(data))]
  ir.pca4 <- prcomp(data_PCA, center = TRUE, scale. = TRUE)
  
  ggbiplot(ir.pca4, choices = components, obs.scale = 1, var.scale = 1, 
    ellipse = F, groups = data[, 2], ellipse.prob = 0.95, circle = F, 
    varname.size = 5, varname.adjust = 3, var.axes = T) + geom_point(aes(shape = groups, 
    colour = groups), size = 3) + geom_hline(yintercept = 0, colour = "black", 
    linetype = "longdash", size = 1) + geom_vline(xintercept = 0, colour = "black", 
    linetype = "longdash", size = 1) + # scale_color_discrete(name = '') +
  theme(legend.direction = "horizontal", legend.position = "top") + stat_ellipse(type = "t", 
    size = 1, geom = "polygon", alpha = 0.2, aes(fill = groups), level = 0.85) + 
    stat_ellipse(type = "t", size = 1, aes(colour = groups), level = 0.85) + 
    guides(groups = FALSE) + theme(text = element_text(size = 15)) + 
    ggtitle("PCA") + theme(plot.title = element_text(hjust = 0.5))
}
