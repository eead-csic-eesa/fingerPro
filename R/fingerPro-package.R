#' A comprehensive package for sediment source unmixing
#' @name fingerPro
#'
#' @description
#' This package quantifies the provenance of sediments in a catchment or study area. Based on a characterization of the sediment sources and the end sediment mixtures, a mixing model algorithm is applied to the sediment mixtures to estimate the relative contribution of each potential source. The package includes several graphs to help users in their data understanding, such as box plots, correlation, PCA, and LDA graphs. In addition, new developments such as the Consensus Ranking (CR), Consistent Tracer Selection (CTS), and Linear Variability Propagation (LVP) methods are included to correctly apply the fingerprinting technique and increase dataset and model understanding. A new method based on Conservative Balance (CB) method has also been included to enable the use of isotopic tracers.
#'
#' @author
#' * Borja Latorre
#' * Ivan Lizaga
#' * Leticia Gaspar
#' * Leticia Palazon
#' * Ana Navas
#' * Maintainer: Erosion, and Soil and Water Evaluation (Research Group) <fingerpro@eead.csic.es>
#'
#' @md
#'
#' @seealso
#' Useful links:
#' * [GitHub repository](https://github.com/eead-csic-eesa/fingerPro)
#'
#' @section Legal Deposits:
#'
#' - FingerPro R. An R package for sediment source fingerprinting (computer
#'   program). Authors: Iván Lizaga, Borja Latorre, Leticia Gaspar, Ana 
#'   María Navas. (EEAD-CSIC). Notarial Act No. 3758 (José Periel Martín), 
#'   18/10/2019. Representative of CSIC: Javier Echave Oria.
#'
#' - FingerPro. Model for environmental mixture analysis (computer 
#'   program). Authors: Leticia Palazón, Borja Latorre, Ana María Navas. 
#'   (EEAD-CSIC). Notarial Act No. 4021 (Pedro Antonio Mateos Salgado), 
#'   21/07/2017. Representative of CSIC: Javier Echave Oria.
#'
#' @examples
#' # Load the 'fingerPro' package to access its functions.
#' library('fingerPro')
#' 
#' ################################################################################
#' ##########################  EXPLORATORY DATA ###################################
#' ################################################################################
#' 
#' # Load the example dataset for a 3-source mixing problem.
#' data <- read.csv(system.file("extdata", "example_geo_3s_raw.csv", package = "fingerPro") )
#' 
#' # Verify the structure and integrity of the loaded dataset.
#' check_database(data)
#' 
#' # Create a box and whisker plot to visualize the distribution of each tracer.
#' # box_plot <- box_plot(data)
#' 
#' # Generate a correlation plot to examine relationships between tracers.
#' # correlation_plot(data)
#' 
#' # Perform and plot Linear Discriminant Analysis (LDA) to visualize group separation.
#' # LDA_plot <- LDA_plot(data)
#' 
#' # Perform and plot Principal Component Analysis (PCA) for dimensionality reduction.
#' # PCA_plot(data)
#' 
#' # Perform a Kruskal-Wallis test (KW) (p-values less than 0.05)
#' # output_KW <- KW_test(data, pvalue=0.05)
#' 
#' # Perform Discriminant Function Analysis (DFA) (confidence level is set to 0.1)
#' # output_DFA <- DFA_test(data, niveau = 0.1)
#' 
#' ################################################################################
#' ####################  FingerPro TRACER SELECTION  ##############################
#' ################################################################################
#' 
#' # Individual Tracer Analysis (ITA) to get descriptive statistics for each tracer.
#' # output_ITA <- individual_tracer_analysis(data)
#' 
#' # Calculate the Conservativeness Index (CI) for each tracer based on ITA results.
#' # output_CI <- CI(output_ITA)
#' 
#' # Generate ternary diagrams to visualize source contribution of tracers.
#' # ternary_diagram(output_ITA, tracers = c(1:9), rows = 3, cols = 3, solution = NA)
#' 
#' # Perform the Range Test (RT) to identify non-conservative tracers.
#' # output_RT <- range_test(data)
#' 
#' # Calculates a Consensus Ranking (CR) score to identify the most reliable tracers.
#' # output_CR <- CR(data, debates = 1000)
#' 
#' # Extract all possible minimal tracer combinations (seeds) for further evaluation.
#' # This helps identify the most discriminant subsets of tracers.
#' # output_CTS_seeds <- CTS_seeds(data, iter = 1000)
#' 
#' # Evaluate the mathematical consistency of a specific tracer combination using the
#' # Consistency Test and Selection (CTS) error.
#' # The user must select a row from the CTS_seeds output based on a criteria (see help CTS_error)
#' # Criteria: positive apportionment (if negative close to zero), high percentage of 
#' # physically feasible solutions, low dispersion)
#' # e.g. select row 1: solution = output_CTS_seeds[1, ]
#' # output_CTS <- CTS_error(data, solution = output_CTS_seeds[1,])
#' 
#' #### OPTIMUM TRACER SELECTION
#' 
#' # Merge the results from CTS, CR, and CI into a single summary data frame.
#' # output_data_summary <- merge(output_CTS, output_CR, by = "tracer")
#' # output_data_summary <- merge(output_data_summary, output_CI, by = "tracer")
#' 
#' # Filter the summary data to select only the most robust tracers.
#' # The criteria are a low CTS error (< 0.05) and a high CR score (> 80).
#' # output_data_TracerSelection <- output_data_summary[output_data_summary$CTS_err < 0.05 &
#' #                                                    output_data_summary$CR_score > 80, ]
#' 
#' ################################################################################
#' ########################  U N M I X I N G  #####################################
#' ################################################################################
#' 
#' # Reload the original data to ensure the analysis starts from the full dataset.
#' # data <- read.csv(system.file("extdata", "example_geo_3s_raw.csv", package = "fingerPro") )
#' 
#' # Select only the tracers identified as optimal in the previous step.
#' # data <- select_tracers(data, output_data_TracerSelection[, 1])
#' 
#' # Run the unmixing model to estimate source apportionment.
#' # output_UNMIX <- unmix(data)
#' 
#' # Plot the unmixing results using a violin plot to visualize source contributions.
#' # plot_results(output_UNMIX)
#'  
#' # Plot the unmixing results using a density plot for an alternative visualization.
#' # plot_results(output_UNMIX, violin = FALSE)
#'

#' @keywords internal
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
NULL
