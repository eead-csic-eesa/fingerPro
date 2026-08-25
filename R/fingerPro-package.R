#' Unmixing Model Framework
#' @name fingerPro
#'
#' @description
#' Quantifies the provenance of sediments by applying a mixing model algorithm to end sediment mixtures based on a comprehensive characterization of the sediment sources. The fingerPro model builds upon the foundational concept of using mass balance linear equations for sediment source quantification by incorporating several distinct technical advancements. It employs an optimization approach to normalize discrepancies in tracer ranges and minimize the objective function. Latin hypercube sampling is used to explore all possible combinations of source contributions (0-100%), mitigating the risk of local minima. Uncertainty in source estimates is quantified through a Monte Carlo routine, and the model includes additional metrics, such as the normalized error of the virtual mixture, to detect mathematical inconsistencies, non-physical solutions, and biases. A new linear variability propagation (LVP) method is also included to address and quantify potential bias in model outcomes, particularly when dealing with dominant or non-contributing sources and high source variability, offering a significant advancement for field studies where direct comparison with theoretical apportionments is not feasible. In addition to the unmixing model, a complete framework for tracer selection is included. Several methods are implemented to evaluate tracer behaviour by considering both source and mixture information. These include the Consistent Tracer Selection (CTS) method to explore all tracer combinations and select the optimal ones improving the robustness and interpretability of the model results. A Conservative Balance (CB) method is also incorporated to enable the use of isotopic tracers. The package also provides several graphical tools to support data exploration and interpretation, including box plots, correlation plots, Linear Discriminant Analysis (LDA) and Principal Component Analysis (PCA).
#'
#' @author
#' * Borja Latorre (Core Team)
#' * Leticia Gaspar (Core Team)
#' * Ivan Lizaga
#' * Leticia Palazon
#' * Ana Navas (Core Team)
#' * Maintainer: Erosion and Evaluation of Soil and Water (EESA research Group) <fingerpro@eead.csic.es>
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
#' # Load the example dataset for a 3-source mixing problem.
#' data <- read.csv(system.file("extdata", "example_geochemical_3s_raw.csv", package = "fingerPro") )
#' #' 

#' @keywords internal
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
NULL
