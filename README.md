[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1402029.svg)](https://doi.org/10.5281/zenodo.1402029)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/fingerPro)](https://cranchecks.info/pkgs/fingerPro)
# fingerPro-model: Sediment source fingerprinting

You can cite this package on your work as:

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020. FingerPro: an R package for tracking the provenance of sediment. Water Resources Management 272, 111020. https://doi.org/10.1016/j.jenvman.2020.111020**.
and also refer to the code:
**Lizaga I., Latorre B., Gaspar L., Navas A., (2018) fingerPro: An R package for sediment source tracing, https://doi.org/10.5281/zenodo.1402029**.

Soil loss by erosion processes is one of the largest challenges for food production and reservoir siltation around the world. Information on sediment, nutrients and pollutants is required for designing effective control strategies. The estimation of sediment sources is difficult to get using conventional techniques, but sediment fingerprinting is a potentially valuable tool. In this regard, FingerPro includes the procedure for the selection of the best combination of tracers and estimate the contribution of each possible source.

Installation
------------
``` r
# From CRAN
install.packages("fingerPro")
library(fingerPro)
# From fingerPro_1.1.zip
setwd("C:/your/file/directory")
install.packages('fingerPro_1.1.zip', repos = NULL)
```
Example Usage
-------------
``` r
#If you want to use your own data
#setwd("the directory that contains your dataset")
#data <- read.table('your dataset', header = T, sep = '\t')

#Example of the data included in the fingerPro package
#Load the dataset called "catchment" 
install.packages("fingerPro")
library(fingerPro)
# "Catchment": This dataset has been selected from a Mediterranean catchment and contains high-quality radionuclides and geochemistry data.
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
data$Land_Use[data$Land_Use == 'PI1'] <- 'PI'
#select the optimum set of tracers by implementing the statistical tests 
data <- rangeTest(data)
data <- KWTest(data)
data <- DFATest(data)
LDAPlot(data, P3D=F)
#Now the optimum tracer properties selected show a good discrimination, so proceed with the unmix function
``` 
![lda comparisson](https://user-images.githubusercontent.com/30837036/43969407-535dab0c-9cca-11e8-8c3d-18fbc3048fb0.png)
```r
# Figure a) LDAPlot previous to PI1+PI fusion. b) LDAPlot after the fusion of both pines
result <- unmix(data, samples = 100L, iter =100L)
#Display the results
plotResults(result, y_high = 5, n = 1)
```

![results](https://user-images.githubusercontent.com/30837036/43969666-2ebd7a1a-9ccb-11e8-8d71-445ad2e15daa.png)
```r
    id   GOF.mean     GOF.SD    AG.mean      AG.SD    PI.mean      PI.SD    SS.mean      SS.SD
  42744 0.94430071 0.03681212  0.18148918  0.061388  0.4726643   0.0785878  0.3458461  0.0654922
    
```
