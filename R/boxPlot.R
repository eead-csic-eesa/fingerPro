#' Box and whiskers plot
#'
#' The boxplot compactly shows the distribution of a continuous variable. It displays five summary statistics (the median, two hinges and two whiskers), and all "outlying" points individually.
#'
#' @param data Data frame containing source and mixtures data
#' @param columns Numeric vector containing the index of the columns in the chart (the first column refers to the first variable)
#' @param ncol Number of charts per row
#'
#' @export
#'
boxPlot <- function(data, columns=1:ncol(data)-2, ncol=3) {
  # reorder the groups
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  dat_plot <- melt(data, id = c(1, 2))

  funclist <- list()
  for (i in 1:length(columns)) {
    eval(parse(text=paste(
    "funclist[[", toString(i), "]] <- function(x) {\n",
    " y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[columns[", toString(i), "]],4]\n",
    " return(data.frame(y = min(y) - 0.05 * (max(y) - min(y)), label = paste0('n=', length(x))))\n",
    "}"
    )))
  }

  funclist2 <- list()
  for (i in 1:length(columns)) {
    eval(parse(text=paste(
    "funclist2[[", toString(i), "]] <- function(x) {\n",
    "  y <- dat_plot[dat_plot$variable == unique(dat_plot$variable)[columns[", toString(i), "]],4]\n",
    "  return(data.frame(y = median(x) + 0.03 * (max(y) - min(y)), label = round(mean(x), 1)))\n",
    "}"
    )))
  }

  g1 <- ggplot(dat_plot[dat_plot$variable == unique(dat_plot$variable)[columns[1]],], aes_string(x = colnames(data)[2], y = "value", color = colnames(data)[2])) + geom_boxplot()
  g_legend<-function(a.gplot) {
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  mylegend <- g_legend(g1)

  blank <- grid.rect(gp=gpar(col="white"))

  glist <- list()
  j <- 1
  for (i in 1:length(columns)) { 
    g1 <- ggplot(dat_plot[dat_plot$variable == unique(dat_plot$variable)[columns[i]],], aes_string(x = colnames(data)[2], y = "value", color = colnames(data)[2])) + geom_boxplot() +
      ggtitle(unique(dat_plot$variable)[columns[i]]) + theme(plot.title = element_text(hjust = 0.5)) + xlab("") + ylab("") + 
      stat_summary(fun.data = funclist[[i]],  geom = "text", fontface="bold", color="white", alpha=0.5,size=3.9) + 
      stat_summary(fun.data = funclist[[i]],  geom = "text", fontface="bold", color="white", alpha=0.5,size=4.0) + 
      stat_summary(fun.data = funclist[[i]],  geom = "text", fontface="bold", color="white", alpha=0.5,size=4.1) + 
      stat_summary(fun.data = funclist[[i]],  geom = "text", fontface="bold", color="white", alpha=0.5,size=4.2) + 
      stat_summary(fun.data = funclist[[i]],  geom = "text") + 
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface="bold", color="white", alpha=0.5,size=3.9) + 
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface="bold", color="white", alpha=0.5,size=4.0) + 
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface="bold", color="white", alpha=0.5,size=4.1) + 
      stat_summary(fun.data = funclist2[[i]], geom = "text", fontface="bold", color="white", alpha=0.5,size=4.2) + 
      stat_summary(fun.data = funclist2[[i]], geom = "text") + 
      theme(legend.position="none")

    glist[[j]] <- g1

    if (i == ncol) {
      j <- j + 1
      glist[[j]] <- mylegend
    }
    else if (i %% ncol == 0) {
      j <- j + 1
      glist[[j]] <- blank
    }

    j <- j + 1
  }

  p3 <- grid.arrange(grobs = glist, ncol=ncol+1, widths=c(rep(10, each=ncol), 3))
}

