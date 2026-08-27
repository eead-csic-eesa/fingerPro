#' @title Generate box-and-whisker plots for sediment tracers
#' @description This function creates a series of box and whisker plots arranged in a grid.
#' It uses a paging system to prevent overlapping and ensures equal-sized plots.
#'
#' @param data A data frame containing sediment source and mixture data.
#' @param page Integer specifying which set of tracers to display (default = 1).
#' @param n_row Number of rows per page (default = 3).
#' @param n_col Number of columns per page (default = 2).
#' @param colors Optional character vector of colors for the groups.
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import reshape
#' @export
box_plot <- function(data, page = 1, n_row = 2, n_col = 3, colors = NULL) {

  # 1. Initial Data Prep (fingerPro standards)
  if(is_averaged(data)) {
    data <- raw_dataset(data)
  }
  
  # Ensure groups are factors for consistent plotting
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))
  group_col <- colnames(data)[2]
  
  # Melt data for ggplot2
  dat_plot <- reshape::melt(data, id = c(1, 2))
  
  # 2. Paging Logic
  all_tracer_names <- unique(dat_plot$variable)
  total_tracers <- length(all_tracer_names)
  plots_per_page <- n_row * n_col
  total_pages <- ceiling(total_tracers / plots_per_page)

  if (page > total_pages || page < 1) {
    stop(paste0("Page ", page, " not found. Total pages available: ", total_pages))
  }
  
  start_idx <- ((page - 1) * plots_per_page) + 1
  end_idx <- min(page * plots_per_page, total_tracers)
  selected_tracers <- all_tracer_names[start_idx:end_idx]

  # Console Feedback
  message(sprintf("Page %i/%i | Showing tracers %i to %i", 
                  page, total_pages, start_idx, end_idx))

  # 3. Legend Extraction Utility
  g_legend <- function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    if (length(leg) > 0) return(tmp$grobs[[leg]]) else return(NULL)
  }

  # 4. Generate the List of Plots
  glist <- list()
  mylegend <- NULL

  for (tracer in selected_tracers) {
    
    p_item <- local({
      
      current_tracer <- tracer
      plot_sub <- dat_plot[dat_plot$variable == current_tracer, ]
      
      # Calculate range and min for THIS specific tracer graph
      y_vals <- plot_sub$value
      y_min_total <- min(y_vals, na.rm = TRUE)
      y_max_total <- max(y_vals, na.rm = TRUE)
      y_range_total <- y_max_total - y_min_total
      
      # Buffer for the labels at the bottom
      plot_bottom <- y_min_total - (0.05 * y_range_total)

      p <- ggplot(plot_sub, aes(x = .data[[group_col]], y = value, color = .data[[group_col]])) +
        geom_boxplot() +
        ggtitle(current_tracer) +
        # Force the y-axis to accommodate the labels at the bottom consistently
        expand_limits(y = plot_bottom) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          panel.grid.minor = element_blank()
        ) +
        xlab("") + ylab("") +
        
        # Add n= sample size at the bottom (y_min_total is equal for all groups here)
        stat_summary(fun.data = function(x) {
          data.frame(y = plot_bottom, label = paste0('n=', length(x)))
        }, geom = "text", size = 3.7, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = plot_bottom, label = paste0('n=', length(x)))
        }, geom = "text", size = 3.4, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = plot_bottom, label = paste0('n=', length(x)))
        }, geom = "text", size = 3.1, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = plot_bottom, label = paste0('n=', length(x)))
        }, geom = "text", size = 2.8) +
        
        # Add Mean value above the median
        stat_summary(fun.data = function(x) {
          data.frame(y = median(x) + 0.05 * y_range_total, label = round(mean(x), 1))
        }, geom = "text", size = 4.1, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = median(x) + 0.05 * y_range_total, label = round(mean(x), 1))
        }, geom = "text", size = 3.8, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = median(x) + 0.05 * y_range_total, label = round(mean(x), 1))
        }, geom = "text", size = 3.5, fontface = "bold", color="white", alpha = 0.5) + 
        stat_summary(fun.data = function(x) {
          data.frame(y = median(x) + 0.05 * y_range_total, label = round(mean(x), 1))
        }, geom = "text", size = 3.2)

      # Apply manual colors if provided
      if (!is.null(colors)) {
        p <- p + scale_color_manual(values = colors)
      }
      
      return(p)
    })
    
    # Capture the legend from the first plot processed
    if (is.null(mylegend)) {
      mylegend <- g_legend(p_item + theme(legend.position = "right", legend.title = element_blank()))
    }
    
    glist[[length(glist) + 1]] <- p_item
  }

  # 4. Final Assembly (No Legend)
  # Simply arrange the grobs in the specified grid without adding a legend column
  gridExtra::grid.arrange(
    grobs = glist, 
    nrow = n_row, 
    ncol = n_col
  )
}
