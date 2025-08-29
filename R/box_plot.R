#' @title Box and whiskers plot for sediment tracers
#'
#' @description This function creates a series of box and whisker plots to visualize the distribution and variability of individual tracers within a dataset. It is designed to work with sediment source and mixture data, automatically adapting to averaged or raw data formats.
#'
#' @param data A data frame containing sediment source and mixture data. Users should ensure their data is in a valid format by using the check_database() function before running this function.
#' @param tracers A numeric vector specifying the column indices of the tracers to be plotted. The index 1 corresponds to the first tracer column after the sample ID and group columns. If NULL (the default), plots will be generated for all tracer columns.
#' @param ncol An integer specifying the number of charts to display per row in the final plot layout.
#' @param colors A character vector of colors to use for the box plots. The colors are applied sequentially to each group (sources and mixture).
#'
#' @details This function is a wrapper for ggplot2 that automates the creation of a series of box plots, one for each tracer. The function first checks if the input data is averaged, and if so, converts it to a virtual raw dataset using the raw_dataset() function to enable the box plot visualization.
#'
#' Each plot displays the distribution of a single tracer, with different groups (sources and mixtures) represented by separate box plots. In addition to the standard five-number summary (median, hinges, and whiskers), the function also overlays the sample count and the mean value for each group, providing a more detailed summary of the data.
#'
#' The final output is a multi-panel plot arranged in a grid, with an optional legend depending on the input data.
#'
#' @export
box_plot <- function(data, tracers = NULL, ncol = 3, colors = NULL) {

	# If data is averaged, convert it to a raw dataset
	if(is_averaged(data)) {
		data <- raw_dataset(data)
	}
	
	if(is.null(tracers)) {
		tracers <- c(1:(ncol(inputMixture(data))-1))
	}

  # reorder the groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  dat_plot <- reshape::melt(data, id = c(1, 2))

  funclist <- list()
  for (i in 1:length(tracers)) {
    eval(parse(text = paste(
      "funclist[[", toString(i), "]] <- function(x) {\n",
      " y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[tracers[", toString(i), "]],4]\n",
      " return(data.frame(y = min(y) - 0.05 * (max(y) - min(y)), label = paste0('n=', length(x))))\n",
      "}"
    )))
  }

  funclist2 <- list()
  for (i in 1:length(tracers)) {
    eval(parse(text = paste(
      "funclist2[[", toString(i), "]] <- function(x) {\n",
      "  y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[tracers[", toString(i), "]],4]\n",
      "  return(data.frame(y = median(x) + 0.03 * (max(y) - min(y)), label = round(mean(x), 1)))\n",
      "}"
    )))
  }

  g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) {
      legend <- tmp$grobs[[leg]]
      return(legend)
    } else {
      return(NULL)
    }
  }

  blank <- grid::grid.rect(gp = grid::gpar(col = "white"))

  glist <- list()
  j <- 1
  for (i in 1:length(tracers)) {
    p <- ggplot(dat_plot[dat_plot$variable == unique(dat_plot$variable)[tracers[i]], ],
                aes_string(x = colnames(data)[2], y = "value", color = colnames(data)[2])) +
      geom_boxplot() +
      ggtitle(unique(dat_plot$variable)[tracers[i]]) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlab("") +
      ylab("") +
      stat_summary(fun.data = funclist[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 3.9) +
      stat_summary(fun.data = funclist[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.0) +
      stat_summary(fun.data = funclist[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.1) +
      stat_summary(fun.data = funclist[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.2) +
      stat_summary(fun.data = funclist[[i]], geom = "text", aes(color = !!sym(colnames(data)[2]))) +
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 3.9) +
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.0) +
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.1) +
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface = "bold", color="white", alpha = 0.5, size = 4.2) +
      stat_summary(fun.data = funclist2[[i]], geom = "text", aes(color = !!sym(colnames(data)[2]))) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme(legend.position = "none")
    
    if (!is.null(colors)) {
      p <- p + scale_color_manual(values = colors) + scale_fill_manual(values = colors)
    }
    
    glist[[j]] <- p
    
    if (i == 1) {
      mylegend <- g_legend(p)
      if (!is.null(mylegend)) {
        j <- j + 1
        glist[[j]] <- mylegend
      }
    } else if (i %% ncol == 0) {
      j <- j + 1
      glist[[j]] <- blank
    }
    
    j <- j + 1
  }
  
  p3 <- gridExtra::grid.arrange(grobs = glist, ncol = ncol + 1, widths = c(rep(10, each = ncol), 3))
}
