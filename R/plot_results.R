#' Displays the results in the screen
#'
#' The function performs a density chart of the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#' 
#' @param data Data frame containing the relative contribution of the potential sediment sources for each sediment mixture in the dataset
#' @param y_high Number of the vertical height of the y-axis
#' @param n Number of charts per row
#' @param scaled Boolean to switch between regular density plot to scaled density plot
#' @param violin Boolean to switch between density to violin chart
#' @param x_limit Number of the horizontal height of the y-axis
#' @param colors Vector of colors to use for the plots
#' 
#' @export
#' 
plotResults <- function(data, y_high = 1, n = 1, scaled = T, violin = F, colors = NULL) {
  data_plots <- melt(data, id = c(1:2))
  
  if (violin == F) {  
    if (scaled == T) {
      plot <- ggplot(data_plots, aes(x = value)) + 
        geom_density(aes(group = variable,  colour = variable, fill = variable, y = after_stat(scaled)), n = 200, bw = 0.055, alpha = 0.35) + 
        # xlim(-0.2, 1.2) + #-0.5,1.5
        coord_cartesian(xlim = c(-0.2, 1.2), expand = T)
    } else {
      plot <- ggplot(data_plots, aes(x = value)) + 
        geom_density(aes(group = variable, colour = variable, fill = variable), n = 200, bw = 0.06, alpha = 0.35) +
        coord_cartesian(xlim = c(0, 1)) + #-0.5,1.5
        ylim(0, y_high)
    }
    
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) + scale_color_manual(values = colors)
    } else {
      plot <- plot + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    }
    
    # Additional lines retained:
    # +
    #   geom_density(data = data_plots1, aes(group = variable, colour = variable ,fill = variable), size= 1 )
    
    # plot_GOF <- ggplot(data_plots, aes(x = GOF)) + geom_density(aes(colour = GOF, fill = "red"), alpha = 0.35) +
    #   xlim(0, 1) +  scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    
    # plot <- ggplot(data_plots, aes(x = value)) +
    # geom_freqpoly(aes(group = variable, colour = variable), binwidth = 0.05) +
    # xlim(0, 1)
    # plot <- ggplot(data_plots, aes(x = value)) +
    # geom_histogram(aes(group = variable, colour = variable, fill = variable),
    # alpha = 0.2, binwidth = 0.025) + xlim(0, 1)
  } else {  
    plot <- ggplot(data_plots, aes(variable, value, colour = variable)) + 
      geom_violin(bw = 0.06, trim = F, alpha = 0.65) +
      geom_boxplot(width = 0.2, alpha = 0.65, colour = "black", aes(fill = variable)) +
      labs(title = "Violin Plot of length by contribution", x = "Source", y = "apportionment (%)") +
      stat_summary(fun = mean, geom = "point", size = 2.5, color = "red", shape = 8) +  
      coord_cartesian(ylim = c(0, 1), expand = T) #-0.5,1.5
    
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) + scale_color_manual(values = colors)
    } else {
      plot <- plot + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    }
  }
  
  plot <- plot + facet_wrap(~id, ncol = n)
  
  # Additional lines retained:
  # dat <- print(head(data, 1))
  # plotgof <- plot_GOF + guides(fill = FALSE)
  # plotgof <- plotgof + facet_wrap(~id, ncol = n)
  # print(plotgof)
  
  print(plot)
}
