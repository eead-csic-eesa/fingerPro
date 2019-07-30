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

Maintainer: Ivan Lizaga <ilizaga@eead.csic.es // lizaga.ivan10@gmail.com>
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{ Sediment source fingerprinting package }
\seealso{
  https://github.com/eead-csic-eesa
}
\examples{
#Created by Ivan Lizaga 30/07/2019

#If you want to use your own data
#setwd("the directory that contains your dataset")
#data <- read.table('your dataset', header = T, sep = '\t')

#Example of the data included in the fingerPro package
#Load the dataset called "catchment" 
install.packages("fingerPro")
library(fingerPro)
# "Catchment": this dataset has been selected from a Mediterranean catchment for this purpose and contains high-quality radionuclides and geochemistry data.
#AG (cropland)
#PI and PI1 (Pine forest, at first looks different but when you display de LDA plot you will see that the wisher decision in join both pines as the same source)
#SS (subsoil)
data <- catchment
boxPlot(data, columns = 1:6, ncol = 3)
correlationPlot(data, columns = 1:5, mixtures = T)
LDAPlot(data, P3D=FALSE)
#variables are collinear
#select the optimum set of tracers by implementing the statistical tests 
data <- rangeTest(data)
data <- KWTest(data)
data <- DFATest(data)
#Check how the selected tracers discriminate between sources
LDAPlot(data, P3D=F)
LDAPlot(data, P3D=T)
#2D and 3D LDAPlots suggest that two of the sources have to be combined
#reload the original dataset "catchment"
data <- catchment
# Combine sources PI1 and PI based on the previous LDAPlot
data <- merge_sources(data, "PI1", "PI")
#select the optimum set of tracers by implementing the statistical tests 
data <- rangeTest(data)
data <- KWTest(data)
data <- DFATest(data)
LDAPlot(data, P3D=F)
PCAPlot(data)
#Now the optimum tracer properties selected discriminate well, so proceed with the unmix function
result <- unmix(data, samples = 100L, iter =100L)
#Display the results
plotResults(result, y_high = 5, n = 1, scaled = F)
# Scaled
plotResults(result, y_high = 1, n = 1, scaled = T)
# Violin plots
plotResults(result, n = 2, violin = T)
writeResults(result)
}
