#' Linear discriminant analysis chart
#'
#' The function performs a linear discriminant analysis and displays the data in the relevant dimensions.
#'
#' @param data Data frame containing source and mixtures data
#' @param P3D Boolean to switch between 2 to 3 dimensional chart
#' @param text Boolean to show or not the identification number of each sample point in the plot
#' @param colors Allows choosing between a different set of colors in the plots
#' @param interactive Boolean to determine whether the plot should be interactive
#' 
#' @export
#' 

# Function for Linear Discriminant Analysis (LDA) plot
LDAPlot <- function(data, P3D = FALSE, text = FALSE, colors = NULL, interactive = FALSE) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  
  # read groups (second column)
  groups <- data[, 2]
  
  # assume the last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]
  
  # read sources
  sources <- data[!groups == mixture, ]
  
  # remove mixture level
  s_groups <- droplevels(sources[, 2])
  
  # extract properties
  data.lda1 <- sources[3:ncol(sources)]
  
  # assign groups
  data.lda1$groups <- as.factor(s_groups)
  
  # Check for NA values filled
  if (any(is.na(data.lda1))) {
    message("NA values were filled with the minimum value of the same group.")
  }
  
  # Fill NA values with minimum value of the same group
  data.lda1 <- data.lda1 %>%
    group_by(groups) %>%
    mutate(across(everything(), ~ ifelse(is.na(.), min(., na.rm = TRUE), .)))
  
  # Perform LDA
  data.lda1$groups <- as.factor(s_groups)
  data.lda <- lda(groups ~ ., data = data.lda1) 
  data.lda.pred <- predict(data.lda) 
  data.lda.temp <- data.frame(data.lda.pred$x, Sources = data.lda1$groups)
  
  # Check if 3D plot is requested
  if (P3D == TRUE) {
    # Check if custom colors are provided
    if (!is.null(colors)) {
      # Create a 3D scatter plot with custom colors and transparency
      if (interactive == TRUE) {
        plot <- plot_ly(data.lda.temp, x = ~LD1, y = ~LD2, z = ~LD3, type = "scatter3d", 
                        mode = "markers", color = ~Sources, colors = colors, text = ~Sources,
                        marker = list(size = 8, opacity = 0.7)) %>%
          add_markers(opacity = 0.7)
      } else {
        plot <- with(data.lda.temp, scatter3d(LD1, LD2, LD3, revolution = 1, col = group_number, 
                                              point.col = group_number, speed = 8, groups = Sources, 
                                              bg.col = "white", model.summary = T, surface.alpha = 0.2, 
                                              ellipsoid = TRUE, ellipsoid.alpha = 0.3, level = 0.8))
      }
    } else {
      # Create a 3D scatter plot with default colors and transparency
      if (interactive == TRUE) {
        plot <- plot_ly(data.lda.temp, x = ~LD1, y = ~LD2, z = ~LD3, type = "scatter3d", 
                mode = "markers", color = ~Sources, colors = colors, text = ~Sources,
                marker = list(size = 8, opacity = 0.7)) %>%
          add_markers() %>% add_ellipsoids(x = ~LD1, y = ~LD2, z = ~LD3,
                         opacity = 0.3,     # Adjust the opacity of the ellipsoids
                         line = list(color = "black"))  # Adjust the line color of the ellipsoids
        
        
        
      } else {
        plot <- with(data.lda.temp, scatter3d(LD1, LD2, LD3, revolution = 1, col = group_number, 
                                              point.col = group_number, speed = 8, groups = Sources, 
                                              bg.col = "white", model.summary = T, surface.alpha = 0.2, 
                                              ellipsoid = TRUE, ellipsoid.alpha = 0.3, level = 0.8))
      }
    }
  } else {
    # If 2D plot is requested
    if (text == TRUE) {
      # Add row numbers to the data frame for text labels
      data.lda.temp$row_num = 1:nrow(data.lda.temp)
      # Create a ggplot with text labels
      plot <- ggplot(data.lda.temp, aes(LD1, LD2, colour = Sources, fill = Sources)) +
        geom_point(aes(shape = Sources), size = 3, alpha = 0.7) + 
        geom_text(aes(label = row_num, hjust = 2, vjust = 0.5), size=4) + 
        geom_hline(yintercept = 0, colour = "black", linetype = "longdash") + 
        geom_vline(xintercept = 0, colour = "black", linetype = "longdash") + 
        stat_ellipse(type = "t", size = 1, geom = "polygon", alpha = 0.2,
                     level = 0.9) + 
        stat_ellipse(type = "t", size = 1, alpha = 0.7, level = 0.9) + 
        ggtitle("LDA") + 
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = colors) +   # Specify the colors
        scale_fill_manual(values = colors)      # Specify the colors
    } else {
      # Create a ggplot without text labels
      plot <- ggplot(data.lda.temp, aes(LD1, LD2, colour = Sources, fill = Sources)) + 
        geom_point(aes(shape = Sources), size = 3, alpha = 0.7) + 
        geom_hline(yintercept = 0,colour = "black", linetype = "longdash") + 
        geom_vline(xintercept = 0, colour = "black", linetype = "longdash") + 
        stat_ellipse(type = "t", size = 1, geom = "polygon", alpha = 0.2,
                     level = 0.9) + 
        stat_ellipse(type = "t", size = 1, alpha = 0.7, level = 0.9) + 
        ggtitle("LDA") + 
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = colors) +   # Specify the colors
        scale_fill_manual(values = colors)      # Specify the colors
    }
  }
  
  # Print the resulting plot
  print(plot)
}
