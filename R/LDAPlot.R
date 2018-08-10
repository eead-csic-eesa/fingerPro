#' Linear discriminat analysis chart
#'
#' The function performs a linear discriminant analysis and displays the data in the relevant dimensions.
#'
#' @param data Data frame containing source and mixtures data
#' @param P3D Boolean to switch between 2 to 3 dimensional chart
#'
#' @export
#' 

if(getRversion() >= "2.15.1")  utils::globalVariables(c("LD1", "LD2","LD3","group_number","Sources"))

LDAPlot <- function(data, P3D = FALSE) {
  # reorder factor levels in order of appearance
  data[, 2] <- factor(data[, 2], levels = unique(data[, 2]))

  # read groups (second column)
  groups <- data[, 2]

  # asume last group is mixtures
  mixture <- levels(groups)[nlevels(groups)]

  # read sources
  sources <- data[!groups == mixture, ]

  # remove mixture level
  s_groups <- droplevels(sources[, 2])

  # extract properties
  data.lda1 <- sources[3:ncol(sources)]
  # assign groups
  data.lda1$groups <- as.factor(s_groups)
  data.lda <- lda(groups ~ ., data = data.lda1) 
  data.lda.pred <- predict(data.lda) 
  data.lda.temp <- data.frame(data.lda.pred$x, Sources = data.lda1$groups)

  if (P3D == TRUE) {
    plot <- with(data.lda.temp, scatter3d(LD1, LD2, LD3, revolution = 1, col = group_number, 
      point.col = group_number, speed = 3, groups = Sources, bg.col = "white", 
      model.summary = T, surface.alpha = 0.2, ellipsoid = TRUE, ellipsoid.alpha = 0.3, 
      level = 0.8))  ####defined the 3D LDA
  } else {
    plot <- ggplot(data.lda.temp, aes(LD1, LD2, colour = Sources)) + 
      geom_point(aes(shape = Sources), size = 3) + geom_hline(yintercept = 0, 
      colour = "black", linetype = "longdash") + geom_vline(xintercept = 0, 
      colour = "black", linetype = "longdash") + stat_ellipse(type = "t", 
      size = 1, geom = "polygon", alpha = 0.2, aes(fill = Sources), 
      level = 0.9) + stat_ellipse(type = "t", size = 1, aes(colour = Sources), 
      level = 0.9) + ggtitle("DFA") + theme(plot.title = element_text(hjust = 0.5))
  }
  {
    print(plot)
  }
}
