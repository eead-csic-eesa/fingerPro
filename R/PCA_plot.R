#' Principal component analysis chart
#'
#' The function performs a principal components analysis on the given data matrix and displays a biplot using vqv.ggbiplot package of the results for each different source to help the user in the decision.
#'
#' @param data Data frame containing source and mixtures data
#' @param components Numeric vector containing the index of the two principal components in the chart
#' @param colors Vector of colors to use for the groups in the plot
#'
#' @export
#' 
PCA_plot <- function(data, components = c(1, 2), colors = NULL) {

	# If data is averaged, convert it to a raw dataset
	if(is_averaged(data)) {
		data <- raw_dataset(data)
	}
	
  # Reorder groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  # Extract data for PCA
  data_PCA <- data[, 3:ncol(data)]
  
  # Perform PCA
  ir.pca4 <- prcomp(data_PCA, center = TRUE, scale. = TRUE)
  
  # Extract PCA results
  pca_results <- as.data.frame(ir.pca4$x[, components])
  pca_results$groups <- data[, 2]
  
  # Check for NA values filled
  if (any(is.na(pca_results))) {
    message("NA values were filled with the minimum value of the same group.")
  }
  
  # Fill NA values with minimum value of the same group
  pca_results <- pca_results %>%
    group_by(groups) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), min(., na.rm = TRUE), .)))
  
  # Plot using ggbiplot
  p <- ggbiplot(ir.pca4, choices = components, obs.scale = 1, var.scale = 1, 
                ellipse = FALSE, groups = pca_results$groups, ellipse.prob = 0.95, circle = FALSE, 
                varname.size = 5, varname.adjust = 3, var.axes = TRUE) +
    geom_point(aes(shape = groups, colour = groups), size = 3) +
    geom_hline(yintercept = 0, colour = "black", linetype = "longdash", size = 1) +
    geom_vline(xintercept = 0, colour = "black", linetype = "longdash", size = 1) +
    theme(legend.direction = "horizontal", legend.position = "top") +
    stat_ellipse(type = "t", size = 1, geom = "polygon", alpha = 0.2, aes(fill = groups), level = 0.85) + 
    stat_ellipse(type = "t", size = 1, aes(colour = groups), level = 0.85) + 
    guides(groups = "none") +
    theme(text = element_text(size = 15)) + 
    ggtitle("PCA") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Apply custom colors if provided
  if (!is.null(colors)) {
    # Check if number of colors matches number of unique groups
    num_groups <- length(unique(pca_results$groups))
    if (length(colors) < num_groups) {
      stop("Insufficient number of colors provided.")
    }
    
    p <- p +
      scale_colour_manual(values = colors) +
      scale_fill_manual(values = colors)
  }
  
  return(p)
}
