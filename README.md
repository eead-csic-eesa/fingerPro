[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1402029.svg)](https://doi.org/10.5281/zenodo.1402029)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/fingerPro)](https://cranchecks.info/pkgs/fingerPro)

![logo def-github-04](https://user-images.githubusercontent.com/30837036/91882995-13c90200-ec84-11ea-9643-0191dfbca995.jpg)


FingerPro is an R package that provides you with the tools to i) characterise the different sediment sources, establish correlations between the tracers and assist the selection of the optimal tracers, ii) graph the results using the state of the art of R packages iii) unmix sediment samples to estimate the apportionment of the sediment sources and iv) test the model using data from a Mediterranean study catchment included in the package. 

For additional details, please see the recently published [FingerPro paper](https://link.springer.com/article/10.1007%2Fs11269-020-02650-0) that includes:

- Full description of functions and how to use them
- Full description of equations

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
The following example has also been stored in the package and can be accessed by using the command **?fingerPro** in your R console after the installation step.
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

![plotresults_readme](https://user-images.githubusercontent.com/30837036/91887212-c43a0480-ec8a-11ea-9315-06c666da5507.png)
```r
     id   GOF.mean     GOF.SD   AG.mean     AG.SD    PI.mean      PI.SD    SS.mean      SS.SD
  42744 0.95232049 0.03197755 0.1947956 0.0636432 0.46714587 0.06078968 0.33805860 0.06227879
    
```

### Feedback

This software has been improved by the questions, suggestions, and bug reports of the user community. If you have any comment, please use the [Issues](https://github.com/eead-csic-eesa/fingerPro/issues) page or report them to lizaga.ivan10@gmail.com.


# Citing FingerPro:
You can cite this package on your work as:

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020. FingerPro: an R package for tracking the provenance of sediment. Water Resources Management 272, 111020. https://doi.org/10.1007/s11269-020-02650-0**.

and also refer to the code as:

**Lizaga I., Latorre B., Gaspar L., Navas A., (2018) fingerPro: An R package for sediment source tracing, https://doi.org/10.5281/zenodo.1402029**.

# Related research:

- New tools for understanding individual tracers and [tracer selection methodologies ](https://www.sciencedirect.com/science/article/pii/S0048969720310482?via%3Dihub)
- Sediment source fingerprinting in [Glacial Landscapes](https://www.sciencedirect.com/science/article/pii/S0169555X20302762?via%3Dihub)
- Agricultural Cycle influence in [sediment and pollutant transport](https://www.sciencedirect.com/science/article/pii/S0301479720309488?via%3Dihub)
- Changes in source contribution [during](https://www.sciencedirect.com/science/article/pii/S0301479719304220?via%3Dihub) an exceptional storm event and [before and after the event](https://www.sciencedirect.com/science/article/pii/S0169555X19302302?via%3Dihub)
- [Testing FingerPro model](https://www.sciencedirect.com/science/article/pii/S0016706118300570) with artificial samples
