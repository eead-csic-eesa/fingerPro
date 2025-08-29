#' Displays the results of an unmixing analysis
#'
#' This function generates a plot showing the relative contribution of sediment sources to each mixture. The output of the \code{unmix} function should be used as input for this function.
#' 
#' @param data A data frame, typically the output from the \code{unmix} function, containing the relative contributions of sediment sources.
#' @param violin A logical value. If \code{TRUE}, violin charts are used instead of density plots.
#' @param bounds A numeric vector of length 2 specifying the lower and upper bounds for the data.
#' @param scaled A logical value. If \code{TRUE}, the density plots are scaled.
#' @param y_high The maximum value for the y-axis.
#' @param colors A character vector of colors to use for the plots.
#' @param ncol The number of plots per row.
#' 
#' @export
plot_results <- function(data, violin = T, bounds = c(0.0, 1.0), scaled = T, y_high = 1, colors = NULL, ncol = 1) {
  data_plots <- reshape::melt(data, id = c(1,ncol(data)))
  
  if (violin == F) {  
    if (scaled == T) {
      plot <- ggplot(data_plots, aes(x = value)) + 
        geom_density(aes(group = variable,  colour = variable, fill = variable, y = after_stat(scaled)), bounds = bounds, n = 256, bw = "nrd", alpha = 0.35) + 
        coord_cartesian(xlim = bounds, expand = T)
    } else {
      plot <- ggplot(data_plots, aes(x = value)) + 
        geom_density(aes(group = variable, colour = variable, fill = variable), bounds = bounds, n = 256, bw = "nrd", alpha = 0.35) +
        coord_cartesian(xlim = bounds, expand = T)
        ylim(0, y_high)
    }
    
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) + scale_color_manual(values = colors)
    } else {
      plot <- plot + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    }
  } else {  
    plot <- ggplot(data_plots, aes(variable, value, colour = variable)) + 
      geom_violin(bw = 0.06, trim = F, alpha = 0.65) +
      geom_boxplot(width = 0.2, alpha = 0.65, colour = "black", aes(fill = variable)) +
      labs(x = "Source", y = "Apportionment (%)") +
      stat_summary(fun = mean, geom = "point", size = 2.5, color = "red", shape = 8) +  
      coord_cartesian(ylim = bounds, expand = T)
    
    if (!is.null(colors)) {
      plot <- plot + scale_fill_manual(values = colors) + scale_color_manual(values = colors)
    } else {
      plot <- plot + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    }
  }
  
  plot <- plot + facet_wrap(~ID, ncol = ncol)

  print(plot)
}
