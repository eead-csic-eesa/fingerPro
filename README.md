[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1402029.svg)](https://doi.org/10.5281/zenodo.1402029)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/fingerPro)](https://cranchecks.info/pkgs/fingerPro)
# fingerPro-model: Sediment source fingerprinting

You can cite this package on your work as:
**Lizaga I., Latorre B., Gaspar L., Navas A., (2018) fingerPro: An R package for sediment source tracing, https://doi.org/10.5281/zenodo.1402029**.

Soil erosion is one of the largest challenges for food production and reservoir siltation around the world. Information on sediment, nutrients and pollutants is required for designing effective control strategies. Sediment source estimates are difficult to obtain using traditional monitoring techniques, but sediment source fingerprinting is a potentially valuable tool. This procedure focuses on developing methods that enable the apportionment of sediment sources to be identified from a composite sample of sediment mixtures.
We developed a new tool to quantify the provenance of sediments in a catchment. For the first time the procedure for selection of the best combination of sediment tracers was included in the tool package. A mixing model algorithm is applied to the sediment samples in order to estimate the relative contribution of each potential source. The operations are compiled in an Open Source R package called fingerPro, which unmixes sediment samples after selecting the optimum set of tracers and providing the percentage contribution of each sediment source. A real example was included in the package to test the model. The results highlight the importance of sediment fingerprinting approaches to quantify the main sediment source which will enable us to better understand catchment sediment delivery to water supply reservoirs and the detrimental effect on the quality of the water and aquatic habitats. Refinement of the sediment source fingerprinting technique requires for consistency and transparency in model handling and tracer selection processes.

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
