#' @title Discriminant Function Analysis (DFA) Test
#'
#' @description Performs a stepwise forward variable selection using the Wilk's Lambda criterion to identify the most discriminant tracers in a dataset.
#'
#' @param data A data frame containing the characteristics of sediment sources and mixtures.
#' @param niveau A numeric value specifying the significance level for the approximate F-test decision.
#'
#' @return A data frame containing only the tracers that pass the DFA test.
#'
#' @export
DFA_test <- function(data, niveau = 0.1) {

	# If data is averaged, convert it to a raw dataset
	if(is_averaged(data)) {
		data <- raw_dataset(data)
	}
	
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
  var_A <- sources[3:ncol(sources)]
  
  # assign groups
  var_A$groups <- as.factor(s_groups)
  
  # perform DFA
  DFA <- greedy.wilks(groups ~ ., data = var_A, niveau = niveau)
  
  # filter selected properties
  var <- as.vector(DFA[[1]][, 1])
  DFA_OFF1 <- data[, 3:ncol(data)]
  DFA_OFF1 <- DFA_OFF1[, !(names(DFA_OFF1) %in% var)]
  var_OFF.DFA <- names(DFA_OFF1)
#  total_DFA <- data[, !(names(data) %in% var_OFF.DFA)]
  total_DFA <- data[, c("ID", "samples", var)] # keep tracer order

  cat(" The variable/variables that remains in your dataset are:", 
      crayon::green(crayon::bold(names(total_DFA[3:ncol(total_DFA)]))), ".\n")
  cat("Attention: ", ncol(DFA_OFF1), "variables were removed from a total of", 
      ncol(data[3:ncol(data)]), ":", crayon::red(crayon::bold(var_OFF.DFA)), ".", "\n")

  #return(list(total_DFA, DFA))
  return(total_DFA)
}

