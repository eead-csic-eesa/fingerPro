![github version](https://img.shields.io/badge/GitHub-1.3-blueviolet?logo=github)
[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1402028.svg)](https://doi.org/10.5281/zenodo.1402029)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/fingerPro)](https://cranchecks.info/pkgs/fingerPro)

![logo def-github-04](https://user-images.githubusercontent.com/30837036/91882995-13c90200-ec84-11ea-9643-0191dfbca995.jpg)
### FingerPro is a collaborative R project aiming to solve and share all modelling concerns around the sediment fingerprinting technique. Join us and discover the insights of your data!!

The package provides users with tools to i) characterise you database ii) assist in the selection of the optimal tracers moving forward from the traditional tracer selection methods (also included), iii) extract the multiple and consensual solutions in your database and iv) unmix sediment samples to estimate the apportionment of the sediment sources.
The package used the state of the art equations and techniques for rigorous unmixing avoiding previous thoughts about only relying on models capacity.

<details>
<summary><strong>Table of Contents</strong></summary>
<!-- toc -->

## Table of Contents
-------------
- [Installation](#installation)
- [Preparing your data](#preparing-your-data)
- [Visual plots](#visual-plots)
   * [Box and whisker plots](#box-and-whisker-plots)
   * [Correlation plot](#correlation-plot)
   * [2D and 3D Linear discrimination plot](#2d-and-3d-linear-discrimination-plot)
   * [Principal component analysis plot](#principal-component-analysis-plot)
- [Novel methods for tracer selection and data understanding](#novel-methods-for-tracer-selection-and-data-understanding)
   * [Conservativeness index (CI) and Ternary plots](#conservativeness-index-and-ternary-plots)
   * [Consensus Ranking approach (CR)](#consensus-ranking-approach)
   * [Consistency based tracer selection](#consistency-based-tracer-selection)
- [The basics of the technique (Video)](#the-basics-of-the-technique)
- [Citing FingerPro and its tools](#citing-fingerpro-and-its-tools)
- [Contributing and feedback](#contributing-and-feedback)
- [Related research](#related-research)

<!-- tocstop -->

</details>

--------------------------------------------------------------------------------

For additional details, please see the recently published [FingerPro paper](https://link.springer.com/article/10.1007%2Fs11269-020-02650-0) and the newly developed tools such as the [Consensus Ranking (CR)](https://www.sciencedirect.com/science/article/pii/S0048969720310482?via%3Dihub) and the newest developed [Consistency based tracer selection (CTS)](https://www.sciencedirect.com/science/article/pii/S0048969721028758) that includes:

- Full description of functions and how to use them
- Full description of equations

## Installation
``` r
# From CRAN (version 1.1)
install.packages("fingerPro")
library(fingerPro)
# From GitHub (version 1.3)
devtools::install_github("eead-csic-eesa/fingerPro", ref="master", force = T)
# From your computer (version 1.3)
setwd("C:/your/file/directory")
install.packages('fingerPro_1.3.zip', repos = NULL)
```

## Preparing your data
 To use your own data is as easy as to follow the format supplied in the example data included in the fingerPro package
 ```r
 print (catchment)
 ``` 
 When using raw data the following structure needs to be followed
  - The first column must be numeric corresponding to the sample ID
  - The second column correspond to the different sources
  ```r
 head (catchment[,c(1:6)])
     id Land_Use  Pbex K40 Bi214 Ra226
1 42665       AG  9.48 494  31.6  32.9
2 42666       AG 19.12 470  32.2  35.1
3 42667       AG 22.62 513  31.1  28.6
4 42694       AG 24.08 587  32.2  32.9
5 42741       AG  3.38 567  28.2  29.4
6 42770       AG  0.32 586  29.9  31.6
```
  - Your mixtures must be located in the last rows with the same name in the 2nd column and different ID (column 1)
 ```r
 tail (catchment[,c(1:6)])
       id   Land_Use  Pbex K40 Bi214 Ra226
17 42705        PI1 51.88 532 25.40 25.90
18 42707        PI1 33.50 530 26.00 24.00
19 42807         SS  0.00 478 26.10 24.30
20 42811         SS  0.00 736 26.20 28.20
21 42812         SS  0.00 607 26.15 26.25
22 42744 mix.sample 24.58 456 26.80 26.90
23 42745 mix.sample 25.58 457 27.80 25.90
 ```

If you only have mean and SD data, follow the format supplied in the Kamish dataset example
```r
sources.file <- system.file("extdata", "Raigani.csv", package = "fingerPro")
data <- read.csv(sources.file)
print(data)
```
Once you have your data in the appropriate format, load it to your global environment
```r
setwd("C:/Users/.....")
data <- read.table("your dataset.csv", header = T, sep = ',')
```
## Visual plots
The following example displays all the basics commands available in the package to display informative graphs.

If you want to use your own data see the previous section [preparing your data](#preparing-your-data)

``` r
#Load the dataset called "catchment" 
data <- catchment
```
#### Box and whisker plots
``` r
boxPlot(data, columns = 1:3, ncol = 3)
```
![image](https://user-images.githubusercontent.com/30837036/137011530-a1148d2a-d30e-49f2-bb21-f55579db1797.png)

#### Correlation plot
``` r
correlationPlot(data, columns = 1:6, mixtures = T)
```
![image](https://user-images.githubusercontent.com/30837036/137011798-6f7f6111-3e7b-4f9d-a95d-d88b83deef59.png)

#### 2D and 3D Linear discrimination plot 
``` r
LDAPlot(data[, c(1:10)], P3D = F, text = T)
LDAPlot(data[, c(1:10)], P3D = T)
```
![LDAs-01](https://user-images.githubusercontent.com/30837036/137012581-eec2d87d-fb4c-44fb-904c-09bd124e27fa.png)

#### Principal component analysis plot
``` r
PCAPlot(data, components = 1:2)
```
![image](https://user-images.githubusercontent.com/30837036/137013578-11105fe4-148b-407c-af61-6e6108aa0d1c.png)

# Novel methods for tracer selection and data understanding
-------------
Moving forward from the traditionally implemented tracer selection methods proven wrong by recent research, this section explains how to implement state of the art approaches to extract individual tracer information and multiple solutions to assist you in this crucial step.

Major benefits:
- Understanding your dataset
- Are there specific relationships among my tracers?
- Does my dataset has multiple solutions?
- What's leading my model to its results?
- Agreement between different models (e.g. FingerPro & MixSIAR)

## Conservativeness index and Ternary plots
``` r
sources.file <- system.file("extdata", "Raigani.csv", package = "fingerPro")
data <- read.csv(sources.file)

#Compute the CI index and the individual tracer solutions
results_CI <- CI_Method(data, points = 2000, Means = T) # Means = F (When using raw data)

# Plot the individual tracers solution from the 8 first tracers
Ternary_diagram(results_CI, tracers = c(1:6), nrow = 1, ncol = 6)
```
![image](https://user-images.githubusercontent.com/30837036/136755984-3b9daf01-8362-417f-b69e-08f7d48e584d.png)
##### Ternary plots of all possible predictions of each tracer

## Consensus Ranking approach 
![image](https://user-images.githubusercontent.com/30837036/136759797-278516bb-8312-4880-b542-7ee4ee912c97.png)

``` r
mixture <- tail(data, n = 1)
var <- grep("^D", colnames(data))
mixture <- mixture[-c(var)]
row.names(mixture) <- NULL
source <- head(data,-1)
sgeo <- source[, -1]
mgeo <- mixture[, -1]

# When using raw data 
# sgeo <- inputSource(data)
# mgeo <- inputSample(data)

crgeo <- cr_3s(source=sgeo, mixture = mgeo, maxiter = 2000, seed = 1234567)
head(crgeo)

#RESULTS
   tracer score
1       P 96.90
2      Ba 96.45
3      Li 96.00
4       K 95.15
.       .   .
.       .   .
.       .   .
31     Mn  0.85
32     Pb  0.40
33     Cu  0.35
34      V  0.05
```

## Consistency based tracer selection
![image](https://user-images.githubusercontent.com/30837036/136759634-d3ed9262-9cb8-4f0f-a9d1-36f7fbcbc60c.png)

```r
# compute pairs/triplets (depending on your source numbers)
pgeo <- pairs(sgeo, mgeo, iter = 2000, seed = 1234567)
head(pgeo)
     id        w1          w2           w3        Dw1        Dw2        Dw3   cons       Dmax
1 Ba Sr 0.6317626  0.09009448  0.278142958 0.05204120 0.03254703 0.04903574 0.9925 0.05204120
2  P Sr 0.5959349  0.20660832  0.197456756 0.05287735 0.04940192 0.04803014 0.9970 0.05287735
3 Sr Th 0.6673848 -0.02575151  0.358366677 0.06000961 0.04818316 0.06068657 0.2840 0.06068657
4 Sr Ti 0.4723578  0.60848947 -0.080847215 0.06102190 0.06848843 0.04681482 0.0455 0.06848843
5 Sr Mg 0.5115884  0.48090866  0.007502897 0.04285354 0.07065800 0.07463573 0.5350 0.07463573
6  K Sr 0.6340332  0.08271014  0.283256639 0.05470178 0.07517778 0.07182738 0.8545 0.07517778

#Explore those pairs (e.g. Ba & Sr)
sol <- pgeo[pgeo$id=="Ba Sr",]
ctsgeo <- cts_3s(source=sgeo, mixture=mgeo, sol=c(sol$w1, sol$w2, sol$w3))
ctsgeo <- ctsgeo %>% right_join(crgeo, by=c("tracer"))
ctsgeo <- ctsgeo[ctsgeo$err<0.025 & ctsgeo$score>80,]
print(ctsgeo)
  tracer          err score
1     Ba 1.110223e-16 96.45
3      K 4.117438e-03 95.15
4     Li 1.541905e-02 96.00
8     Sr 1.998401e-15 93.80

data1 <- select(data,"id","sources","Ba", "Li", "K", "Sr", "DBa", "DLi", "DK", "DSr", "n")

result_FP_1 <- unmix(data1, samples = 200, iter = 200, Means = T)
P1 <- plotResults(result_FP_1,y_high = 1)

library(MixSIAR)
result_MixSIAR <- MixSIAR_F_SD(data1, chain = "very long", Means = T)
result_TM_MxS <- as.data.frame(result_MixSIAR)
result_TM_MxS <- cbind(id=1,id2 = 1:nrow(result_TM_MxS), result_TM_MxS)
P2 <- plotResults(result_TM_MxS, y_high = 1)
multiplot(P1, P2, cols = 2)
```
![image](https://user-images.githubusercontent.com/30837036/136755551-b3724283-668f-48e5-a803-dbcdd3a174bc.png)
##### FingerPro and MixSIAR models results from one of the consistent solutions extracted from the CTS method

## The basics of the technique
[![Alt text](https://img.youtube.com/vi/LcrM_vLOa_I/0.jpg)](https://www.youtube.com/watch?v=LcrM_vLOa_I)


## Citing FingerPro and its tools
You can cite this package and the new developed tools on your work as:

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020. FingerPro: an R package for tracking the provenance of sediment. Water Resources Management 272, 111020. https://doi.org/10.1007/s11269-020-02650-0**.

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020a. Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. Science of The Total Environment 720, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537**

**Latorre, B., Lizaga, I., Gaspar, L., Navas, A., 2021. A novel method for analysing consistency and unravelling multiple solutions in sediment fingerprinting. Science of The Total Environment 789, 147804. https://doi.org/10.1016/j.scitotenv.2021.147804**

and also refer to the code as:

**Lizaga I., Latorre B., Gaspar L., Navas A., (2018) fingerPro: An R package for sediment source tracing, https://doi.org/10.5281/zenodo.1402029**.

## Contributing and feedback
This software has been improved by the questions, suggestions, and bug reports of the user community. If you have any comment, please use the [Issues](https://github.com/eead-csic-eesa/fingerPro/issues) page or report them to lizaga.ivan10@gmail.com.

## Related research

- New tools for understanding individual tracers and [tracer selection methodologies ](https://www.sciencedirect.com/science/article/pii/S0048969720310482?via%3Dihub)
- Sediment source fingerprinting in [Glacial Landscapes](https://www.sciencedirect.com/science/article/pii/S0169555X20302762?via%3Dihub)
- Agricultural Cycle influence in [sediment and pollutant transport](https://www.sciencedirect.com/science/article/pii/S0301479720309488?via%3Dihub)
- Changes in source contribution [during](https://www.sciencedirect.com/science/article/pii/S0301479719304220?via%3Dihub) an exceptional storm event and [before and after the event](https://www.sciencedirect.com/science/article/pii/S0169555X19302302?via%3Dihub)
- [Testing FingerPro model](https://www.sciencedirect.com/science/article/pii/S0016706118300570) with artificial samples
