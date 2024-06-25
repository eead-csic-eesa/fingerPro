![github version](https://img.shields.io/badge/GitHub-1.3-blueviolet?logo=github)
[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1402028.svg)](https://doi.org/10.5281/zenodo.1402029)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/fingerPro)](https://cranchecks.info/pkgs/fingerPro)

![logo def-github-04](https://user-images.githubusercontent.com/30837036/91882995-13c90200-ec84-11ea-9643-0191dfbca995.jpg)
### FingerPro is a collaborative R project aiming to solve and share all modelling concerns around the sediment fingerprinting technique. Join us and discover the insights of your data!!

The package provides users with tools to i) characterise you database ii) assist in the selection of the optimal tracers moving forward from the traditional tracer selection methods (also included), iii) extract the multiple and consensual solutions in your database and iv) unmix sediment samples to estimate the apportionment of the sediment sources.
The package used state-of-the-art equations and techniques for rigorous unmixing, avoiding previous thoughts about only relying on model capacity.

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
- [Video Tutorials :clapper:](#video-tutorials)
   * [The basics of the technique](#the-basics-of-the-technique)
   * [The FingerPro package: A tutorial with code examples](#the-fingerpro-package)
- [Citing FingerPro and its tools](#citing-fingerpro-and-its-tools)
- [Contributing and feedback](#contributing-and-feedback)
- [Related research](#related-research)

<!-- tocstop -->

</details>

--------------------------------------------------------------------------------

For additional details, please see the recently published [FingerPro paper](https://link.springer.com/article/10.1007%2Fs11269-020-02650-0) and the newly developed tools such as the [Consensus Ranking (CR)](https://www.sciencedirect.com/science/article/pii/S0048969720310482?via%3Dihub) and the newest developed [Consistency based tracer selection (CTS)](https://www.sciencedirect.com/science/article/pii/S0048969721028758) that includes:

- Full description of functions and how to use them
- Full description of equations

If you're working with **stable isotopes** or want to combine them with elemental tracers, check the latest published method [conservative balance (CB)](https://www.sciencedirect.com/science/article/pii/S0048969722019271?via%3Dihub) and its applicability in a [real case study](https://www.sciencedirect.com/science/article/pii/S0022169424001628?via%3Dihub).

--------------------------------------------------------------------------------
##  :point_right: Frequently Asked Questions :arrow_right: [Check them!](FaQ.md)
--------------------------------------------------------------------------------

## Installation
``` r
# From CRAN (version 1.1)
install.packages("fingerPro")
library(fingerPro)

# From your computer (version 1.3)
Download the fingerPro_1.3.zip file from GitHub on your computer
setwd("C:/your/file/directory")
install.packages('fingerPro_1.3.zip', repos = NULL)

# From GitHub (version 1.3)
devtools::install_github("eead-csic-eesa/fingerPro", ref = "master", force = T)
```

## Preparing your data
 To use your own data is as easy as to follow the format supplied in the example data included in the fingerPro package
 ```r
 print (catchment)
 ``` 
 When using raw data, the following structure needs to be followed
  - The first column must be numeric corresponding to the sample ID
  - The second column corresponds to the different sources
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
       id   Land_Use    Pbex    K40   Bi214   Ra226
18  42707        PI1 33.5000 530.00 26.0000 24.0000
19  42807         SS  0.0000 478.00 26.1000 24.3000
20  42811         SS  0.0000 736.00 26.2000 28.2000
21  42812         SS  0.0000 607.00 26.1500 26.2500
22 428112         SS  0.0000 613.07 26.4115 26.5125
23  42745 mix.sample 24.8258 460.56 27.0680 27.1690
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
The following example displays all the basic commands available in the package to display informative graphs.

If you want to use your own data, see the previous section [preparing your data](#preparing-your-data)

``` r
#Load the dataset called "catchment" 
data <- catchment
```
#### Box and whisker plots
``` r
boxPlot(data, tracers = 1:3, ncol = 3 ,colors = c("#993300", "#33CC00", "#336600", "#9933CC", "#0000FF"))
```
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/02ace1d4-0262-4a14-b48a-672548a6f054" width="750">

#### Correlation plot
``` r
correlationPlot(data, columns = 1:6, mixtures = T, colors = c("#993300", "#33CC00", "#336600", "#9933CC", "#0000FF"))
```
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/1af76282-8313-4404-ae88-4ff964f39acd" width="650">

#### 2D and 3D Linear discrimination plot 
``` r
LDAPlot(data[, c(1:10)],  text = T, P3D = F, interactive = F, colors = c("#993300", "#33CC00", "#336600", "#9933CC"))
LDAPlot(data[, c(1:10)],  text = T, P3D = T, interactive = F, colors = c("#993300", "#33CC00", "#336600", "#9933CC"))
LDAPlot(data[, c(1:10)],  text = T, P3D = T, interactive = T, colors = c("#993300", "#33CC00", "#336600", "#9933CC"))
```
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/f8412a20-1359-40db-8d13-e286b4b41695" width="850">

#### Principal component analysis plot
``` r
PCAPlot(data, components = c(1:2), colors = c("#993300", "#33CC00", "#336600", "#9933CC", "#0000FF"))
```
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/efe5cc89-93ae-40ea-86e9-d45797189f8e" width="350">

# Novel methods for tracer selection and data understanding
-------------
Moving forward from the traditionally implemented tracer selection methods proven wrong by recent research, this section explains how to implement state-of-the-art approaches to extract individual tracer information and multiple solutions to assist you in this crucial step.

Major benefits:
- Understanding your dataset
- Are there specific relationships among my tracers?
- Does my dataset have multiple solutions?
- What's leading my model to its results?
- Agreement between different models (e.g. FingerPro & MixSIAR)

## Conservativeness index and Ternary plots
``` r
sources.file <- system.file("extdata", "Raigani.csv", package = "fingerPro")
data <- read.csv(sources.file)

#Compute the CI index and the individual tracer solutions
results_CI <- CI_Method(data, points = 3000, Means = T) # Means = F (When using raw data)

# Plot the individual tracers solution from the eight first tracers
Ternary_diagram(results_CI, tracers = c(1:6), n_row = 1, n_col = 6)
```
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/b4d649e9-d1af-4249-a2be-ea482229c1a5" width="750">

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

crgeo <- cr_ns(source=sgeo, mixture = mgeo, maxiter = 1000, seed = 1234567)
head(crgeo)

#RESULTS
	tracer score
1      P  96.8
2     Li  96.1
3     Ba  96.0
4      K  95.6
.       .   .
.       .   .
.       .   .
31     Cu   0.7
32     Zn   0.7
33     Pb   0.6
34      V   0.1
```

## Consistency based tracer selection
![image](https://user-images.githubusercontent.com/30837036/136759634-d3ed9262-9cb8-4f0f-a9d1-36f7fbcbc60c.png)

```r
# compute pairs/triplets (depending on your source numbers)
pgeo <- pairs(sgeo, mgeo, iter = 1000, seed = 1234567)
head(pgeo)
	id        w1          w2           w3        Dw1        Dw2        Dw3  cons       Dmax
1  P Sr 0.5959349  0.20660832  0.197456756 0.05212140 0.05042700 0.04688297 1.000 0.05212140
2 Ba Sr 0.6317626  0.09009448  0.278142958 0.05230141 0.03341239 0.04975232 0.995 0.05230141
3 Sr Th 0.6673848 -0.02575151  0.358366677 0.05470869 0.04601870 0.05972316 0.275 0.05972316
4 Sr Ti 0.4723578  0.60848947 -0.080847215 0.05934662 0.06579869 0.04640424 0.038 0.06579869
5 Sr Mg 0.5115884  0.48090866  0.007502897 0.04083850 0.07169679 0.07250162 0.540 0.07250162
6  K Sr 0.6340332  0.08271014  0.283256639 0.05456926 0.07371975 0.07234518 0.851 0.07371975
```
##### Thanks to this method, we can see the presence of multiple solutions that would otherwise go undetected. Now is the time to use all the accumulated expertise from statistical methods, lab work, and fieldwork to decide.

```r
# Explore those pairs and multiple solutions (e.g. P-Sr and Ba-Sr)
# Solution 1
sol_1 <- pgeo[pgeo$id=="P Sr",]
ctsgeo_1 <- cts_3s(source = sgeo, mixture = mgeo, sol = c(sol_1$w1, sol_1$w2, sol_1$w3))
ctsgeo_1 <- ctsgeo_1 %>% right_join(crgeo, by = c("tracer"))
ctsgeo_1 <- ctsgeo_1[ctsgeo_1$err < 0.025 & ctsgeo_1$score > 80,]

# Solution 2
sol_2 <- pgeo[pgeo$id=="Ba Sr",]
ctsgeo_2 <- cts_3s(source = sgeo, mixture = mgeo, sol = c(sol_2$w1, sol_2$w2, sol_2$w3))
ctsgeo_2 <- ctsgeo_2 %>% right_join(crgeo, by = c("tracer"))
ctsgeo_2 <- ctsgeo_2[ctsgeo_2$err < 0.025 & ctsgeo_2$score > 80,]

print(cbind(ctsgeo_1, ctsgeo_2))
   tracer          err score 			tracer          err score
      Li 2.176625e-02  96.1    		 	 Ba 1.110223e-16  96.0
       P 2.220446e-16  96.8      		  K 4.117438e-03  95.6
      Sr 3.330669e-16  93.6     		 Li 1.541905e-02  96.1
      Ca 1.337685e-02  93.0     		 Sr 1.998401e-15  93.6

data_sol_1 <- select(data, "id", "sources", "Ba", "Li", "K", "Sr", "DBa", "DLi", "DK", "DSr", "n")
data_sol_2 <- select(data, "id", "sources", "Li", "P", "Sr", "Ca", "DLi", "DP", "DSr", "DCa", "n")

# Let's unmix the multiple solutions
result_FP_1 <- unmix(data_sol_1, samples = 200, iter = 200, Means = T)
result_FP_2 <- unmix(data_sol_2, samples = 200, iter = 200, Means = T)

P_FP_1 <- plotResults(result_FP_1, y_high = 1, colors = c("#CC0000",  "#33CCFF", "#9933CC"))
P_FP_2 <- plotResults(result_FP_2, y_high = 1, colors = c("#CC0000",  "#33CCFF", "#9933CC"))

grid.arrange(P_FP_1, P_FP_2, ncol=2)
```
#### FingerPro model results from the first two consistent solutions extracted from the CTS method
<img src= "https://github.com/brianstock/MixSIAR/assets/30837036/8e1dbfab-6cfa-4984-9fc2-cac9831539b5" width="750">

# :clapper: Video Tutorials

## The basics of the technique
[![Alt text](https://img.youtube.com/vi/LcrM_vLOa_I/0.jpg)](https://www.youtube.com/watch?v=LcrM_vLOa_I)

**Also available on [bilibili](https://www.bilibili.com/video/BV1xZ4y147Ud/?spm_id_from=333.999.0.0)**

## The FingerPro package
[![Alt text](https://img.youtube.com/vi/7HwGcRSO2O8/0.jpg)](https://www.youtube.com/watch?v=7HwGcRSO2O8&ab_channel=UnmixingScience)

**Also available on [bilibili](https://www.bilibili.com/video/BV1mY411M7AW/?spm_id_from=333.999.0.0&vd_source=be6d6fbd57395fc4b904aa65d8bf97d1)**

## Citing FingerPro and its tools
You can cite this package and the new developed tools on your work as:

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020. FingerPro: an R package for tracking the provenance of sediment. Water Resources Management 272, 111020. https://doi.org/10.1007/s11269-020-02650-0**.

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020a. Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. Science of The Total Environment 720, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537**

**Latorre, B., Lizaga, I., Gaspar, L., Navas, A., 2021. A novel method for analysing consistency and unravelling multiple solutions in sediment fingerprinting. Science of The Total Environment 789, 147804. https://doi.org/10.1016/j.scitotenv.2021.147804**

**Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2022. Combined use of geochemistry and compound-specific stable isotopes for sediment fingerprinting and tracing. Science of The Total Environment 832, 154834. https://doi.org/10.1016/j.scitotenv.2022.154834**

and also refer to the code as:

**Lizaga I., Latorre B., Gaspar L., Navas A., (2018) fingerPro: An R package for sediment source tracing, https://doi.org/10.5281/zenodo.1402029**.

## Contributing and feedback
This software has been improved by the questions, suggestions, and bug reports of the user community. If you have any comments, please use the [Issues](https://github.com/eead-csic-eesa/fingerPro/issues) page or report them to lizaga.ivan10@gmail.com.

## Related research
- Combining geochemistry and [isotopic tracers](https://www.sciencedirect.com/science/article/pii/S0048969722019271?via%3Dihub#f0040)
- New tools for understanding individual tracers and [tracer selection methodologies ](https://www.sciencedirect.com/science/article/pii/S0048969720310482?via%3Dihub)
- Sediment source fingerprinting in Glacial Landscapes, [Svalbard](https://www.sciencedirect.com/science/article/pii/S0169555X20302762?via%3Dihub), [Peruvian Andes](https://onlinelibrary.wiley.com/doi/full/10.1002/hyp.14662) and [Antarctica](https://www.sciencedirect.com/science/article/pii/S0169555X24002629?via%3Dihub)
- Agricultural Cycle influence in [sediment and pollutant transport](https://www.sciencedirect.com/science/article/pii/S0301479720309488?via%3Dihub)
- Changes in source contribution [during](https://www.sciencedirect.com/science/article/pii/S0301479719304220?via%3Dihub) an exceptional storm event and [before and after the event](https://www.sciencedirect.com/science/article/pii/S0169555X19302302?via%3Dihub)
- Sediment source fingerprinting in [desert environments](https://www.sciencedirect.com/science/article/pii/S0341816223001029?via%3Dihub)
- Sediment source fingerprinting to track pollutants in [mountainous fluvial environment](https://www.sciencedirect.com/science/article/pii/S004896972303053X?via%3Dihub), [mining areas](https://www.mdpi.com/2071-1050/15/15/11856) and [agroecosystems](https://www.sciencedirect.com/science/article/pii/S0301479720309488?via%3Dihub) 
- [Testing FingerPro model](https://www.sciencedirect.com/science/article/pii/S0016706118300570) with artificial samples
- [Particle size effect](https://www.sciencedirect.com/science/article/pii/S0169555X2200071X?via%3Dihub)
