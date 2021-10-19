#' Displays the results in the screen
#'
#' The function performs a density chart of the relative contribution of the potential sediment sources for each sediment mixture in the dataset.
#' 
#' @param data Data frame containing the relative contribution of the potential sediment sources for each sediment mixture in the dataset
#' @param y_high Number of the vertical height of the y-axis
#' @param n Number of charts per row
#' @param scaled Boolean to switch between regular density plot to scaled density plot
#' @param violin Boolean to switch between density to violin chart
#' 
#' @export
#' 

if(getRversion() >= "2.15.1")  utils::globalVariables(c("data", "value", "variable"))

plotResults <- function(data, y_high = 1, n = 1, scaled = T, violin = F) {
  data_plots <- melt(data, id = c(1:2))
  
  if (violin == F) {  
    if (scaled == T) {
      plot<-ggplot(data_plots, aes(x = value)) + geom_density(aes(group = variable, 
                                                                  colour = variable, fill = variable, y=..scaled..),
                                                              n=200,bw =0.035, alpha = 0.35) + xlim(0, 1) + 
        ylim(0, 1) + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    } else{
      plot <- ggplot(data_plots, aes(x = value)) + geom_density(aes(group = variable, 
                                                                    colour = variable, fill = variable),n=200,bw=0.06,alpha = 0.35) + xlim(0, 1) + 
        ylim(0, y_high) + scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    }
    
    # +
    #   geom_density(data = data_plots1, aes(group = variable, colour = variable ,fill = variable), size= 1 )
    
    # plot_GOF <- ggplot(data_plots, aes(x = GOF)) + geom_density(aes(colour = GOF, fill = "red"), alpha = 0.35) +
    #   xlim(0, 1) +  scale_fill_discrete(name = "sources") + scale_colour_discrete(name = "sources")
    
    # plot <- ggplot(data_plots, aes(x=value)) +
    # geom_freqpoly(aes(group=variable, colour=variable),binwidth = 0.05)+
    # xlim(0, 1)
    # plot <- ggplot(data_plots, aes(x=value)) +
    # geom_histogram(aes(group=variable, colour=variable, fill=variable),
    # alpha=0.2, binwidth = 0.025)+ xlim(0, 1)
  }
  else {  
    plot <- ggplot(data_plots, aes(variable, value, colour=variable)) + 
      geom_violin(bw=0.06, trim = F, alpha = 0.65)+
      geom_boxplot(width=0.2,alpha=0.65, colour="black", aes(fill=variable))+
      labs(title="Violin Plot of length  by contribution",x="Source", y = "apportionment (%)")+
      stat_summary(fun.y=mean, geom="point", size=2.5, color="red", shape=8)+  ylim(0, 1) 
  }
  plot <- plot + facet_wrap(~id, ncol = n)
  # dat <- print(head(data, 1))
  # plotgof <- plot_GOF + guides(fill=FALSE)
  # plotgof <- plotgof  + facet_wrap(~id, ncol = n)
  # print(plotgof)
  print(plot)
  
}
