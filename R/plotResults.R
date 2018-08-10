#' Displays the results in the screen
#'
#' The function performs a density chart of the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#' 
#' @param data Data frame containing the relative contribution of the potential sediment sources for each sediment mixture in the dataset
#' @param y_high Number of the vertical height of the y-axis
#' @param n Number of charts per row
#' 
#' @export
#' 

if(getRversion() >= "2.15.1")  utils::globalVariables(c("result", "value", "variable"))

plotResults <- function(data, y_high = 6.5, n = 1) {
  data_plots <- melt(result, id = c(1:2))
  # data_plots1 <- melt(head(result[1,]), id = c(1:2))
  id_unicos <- length(unique(data_plots$id))
  
plot<-ggplot(data_plots, aes(x = value)) + geom_density(aes(group = variable, 
    colour = variable, fill = variable),n=200,bw=0.06,alpha = 0.35) + xlim(0, 1) + 
    ylim(0, y_high) + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
  # +
  #   geom_density(data = data_plots1, aes(group = variable, colour = variable ,fill = variable), size= 1 )

  # plot_GOF <- ggplot(data_plots, aes(x = GOF)) + geom_density(aes(colour = GOF, fill = "red"), alpha = 0.35) +
  #   xlim(0, 1) +  scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")

  # plot<- ggplot(data_plots, aes(x=value)) +
  # geom_freqpoly(aes(group=variable, colour=variable),binwidth = 0.05)+
  # xlim(0, 1)
  # plot<- ggplot(data_plots, aes(x=value)) +
  # geom_histogram(aes(group=variable, colour=variable, fill=variable),
  # alpha=0.2, binwidth = 0.025)+ xlim(0, 1)

  plot <- plot + facet_wrap(~id, ncol = n)
  # dat <- print(head(result, 1))
  # plotgof <- plot_GOF + guides(fill=FALSE)
  # plotgof <- plotgof  + facet_wrap(~id, ncol = n)
  # print(plotgof)
  print(plot)
  
 }

