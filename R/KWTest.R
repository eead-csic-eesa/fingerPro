#' Kruskal-Wallis rank sum test
#'
#' This function excludes from the original data frame the properties which do not show significant differences between sources.
#'
#' @param data Data frame containing source and mixtures
#' @param pvalue p-value threshold 
#'
#' @return Data frame only containing the variables that pass the Kruskal-Wallis test
#'
#' @export
#'
KWTest <- function(data, pvalue = 0.05) {
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
  var_s <- sources[, 3:ncol(sources)]
  var <- colnames(var_s)
  
  # perform KW test
  KW_p <- c()
  for (i in 1:ncol(var_s)) {
    KW_p <- append(KW_p, kruskal.test(var_s[, i], s_groups)$p.value)
  }
  KW <- data.frame(var, KW_p)
  
  # remove properties under pvalue
  KWOFF <- subset(KW, KW_p > pvalue)
  varOFFKW <- as.vector(KWOFF$var)
  total_KW <- data[, !(names(data) %in% varOFFKW)]
  
  cat("Attention->", nrow(KWOFF), "variables from a total of", ncol(data[3:ncol(data)]), 
    "were removed:", crayon::red(crayon::bold(varOFFKW)), ".", "\n")
  cat(" The variable/variables that remains in your dataset is/are:", 
      crayon::green(crayon::bold(names(total_KW[3:ncol(total_KW)]))), ".")
  
  return(total_KW)
}
