#' Save the results
#'
#' The function saves the results in the workspace file for all the sediment mixture samples and for each sediment mixture sample separately
#'
#'@param data Data frame containing the relative contribution of the potential sediment sources for each sediment mixture in the dataset
#'
#' @export
#'
writeResults <- function(data) {
#summary <- aggregate(. ~ id, data = data, function(x) c(mean = mean(x), SD = sd(x)))
#Save your data in a .csv file
   # write.csv(summary, row.names = F, file = "Results Summary.csv")
  {
    samples <- split(data, data[,1])
    list2env(samples, envir=.GlobalEnv)
#Save your data in a .csv file
    # for (n in names(samples)) {
    #   write.csv(samples[[n]], row.names = F, file = paste(samples[[n]][, 
    #     "id"][1], ".csv", sep = ""))
  }
}

 