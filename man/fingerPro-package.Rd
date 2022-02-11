\name{fingerPro-package}
\alias{fingerPro-package}
\alias{fingerPro}
\docType{package}
\title{
  Sediment Source Fingerprinting
}
\description{
Soil erosion is one of the biggest challenges for food production and reservoirs siltation around the world. Information on sediment, nutrients and pollutant transport is required for effective control strategies. Source estimates are difficult to obtain using traditional monitoring techniques, but sediment source fingerprinting, has been proved to be a valuable tool. Sediment source fingerprinting offers the potential to assess sediment provenance as a basis to develop management plans and prevent erosion.  The procedure focuses on developing methods that enable the apportionment of sediment sources to be identified from a composite sample of sediment mixture material. We developed an R-package as a tool to quantify the provenance of the sediments in a catchment. A mixing model algorithm is applied to the sediment mixture samples in order to estimate the relative contribution of each potential source. The package consists of a set of functions used to:  i) characterise and pre-process the data, select the optimum subset of tracers; ii) unmix sediment samples and quantify the apportionment of each source; iii) assess the effect of the source variability; and iv) visualize and export the results.
}
\author{
Ivan Lizaga, Borja Latorre, Leticia Gaspar, Ana Navas

Maintainer: Ivan Lizaga <lizaga.ivan10@gmail.com>
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ Sediment source fingerprinting package }
\seealso{
  GitHub repository:
  https://github.com/eead-csic-eesa/fingerPro
  
  VIDEO TUTORIALS:
  https://www.youtube.com/watch?v=LcrM_vLOa_I&ab_channel=UnmixingScience
}
\examples{
# Created by Ivan Lizaga 28/01/2022

# If you want to use your own data,
# I recommend following the format of the datasets included in the package
# setwd("the directory that contains your dataset")
# data <- read.table('your dataset.csv', header = T, sep = ',')

# Example of the data included in the fingerPro package
# Load the dataset called "catchment" 
install.packages("fingerPro")
library(fingerPro)
# "catchment": this dataset has been selected from a Mediterranean catchment
# to test the multiple functions inside the package. 
# It contains high-quality radionuclides and geochemistry data, four sources and one mixture.
# AG (cropland), PI and PI1 (Two types of Pine forest), and SS (subsoil)


# Load the 'catchment' dataset to test the graphical functions
data <- catchment

# Display boxplots and a correlation graph of the loaded dataset
boxPlot(data, columns = 1:6, ncol = 3)
correlationPlot(data, columns = 1:7, mixtures = TRUE)

# Check how the selected tracers discriminate between the potential sources
LDAPlot(data[, c(1:10)], P3D=FALSE, text = FALSE)
LDAPlot(data[, c(1:10)], P3D=FALSE, text = TRUE) #adding point information
LDAPlot(data[, c(1:10)], P3D=TRUE, text = FALSE) #3D LDA

# Displaying a PCA graph 
PCAPlot(data)

# Imagine we would like to combine PI and PI1 sources based on the previous graphs
data$Land_Use[data$Land_Use == 'PI1'] <- 'PI'

# Let's check again the discriminant capacity of the first 10 tracers
LDAPlot(data[, c(1:10)], P3D=FALSE, text = TRUE)
}