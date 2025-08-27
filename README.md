![github version](https://img.shields.io/badge/GitHub-2.0-blueviolet?logo=github)
[![cran version](http://www.r-pkg.org/badges/version/fingerPro?color=yellow)](https://cran.r-project.org/package=fingerPro)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/fingerPro)](https://github.com/metacran/cranlogs.app)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

![logo def-github-04](https://user-images.githubusercontent.com/30837036/91882995-13c90200-ec84-11ea-9643-0191dfbca995.jpg)
# A Comprehensive R Package for Sediment Source Unmixing

`fingerPro` is an R package designed to quantify the provenance of sediments in a catchment or study area. By characterizing sediment sources and end sediment mixtures, the package applies a mixing model algorithm to estimate the relative contribution of each potential source. It includes various functions for data understanding, such as box plots, correlation plots, PCA, and LDA graphs. The package also incorporates advanced methods like Consensus Ranking (CR), Consistent Tracer Selection (CTS), and Linear Variability Propagation (LVP) to enhance the application of the fingerprinting technique. It also includes a new Conservative Balance (CB) method to enable the use of isotopic tracers.

-----

### Key Features

#### Data Exploration

  * **Box and whiskers plot**: This function creates a series of plots to visualize the distribution and variability of individual tracers. It works with both averaged and raw data formats.
  * **Correlation matrix chart**: This function displays a correlation matrix of each property, divided by the different sources, to help with tracer selection.
  * **Linear discriminant analysis chart**: This function performs a linear discriminant analysis and displays the data in the relevant dimensions.
  * **Principal component analysis chart**: This function performs a principal components analysis on the data and displays a biplot of the results for each source.

#### Tracer Selection Methods

  * **Consensus Ranking (CR)**: This is an ensemble technique to identify non-conservative and dissenting tracers. The CR score, which ranges from 100 to 0, indicates a tracer's rank in terms of consensus and conservativeness. Tracers are ordered by score, with high scores indicating conservative tracers and low scores indicating dissenting ones.
  * **Conservativeness Index (CI)**: This function calculates the CI for each tracer based on the results of an individual tracer analysis. The CI was adapted from its original definition to better describe tracer conservativeness in a high-dimensional space of multiple sources.
  * **Consistent Tracer Selection (CTS) seeds**: This function extracts all possible minimal tracer combinations to identify the most discriminant ones. The dispersion of the solution reflects the discriminant capacity of each tracer combination, where a lower dispersion indicates higher discriminant capacity.
  * **Range test**: This function excludes properties of the sediment mixture that are outside the minimum and maximum values of the sediment sources.

#### Specialized Tracer Analysis

  * **Conservative Balance (CB) Method**: This function transforms isotopic ratio and content data into virtual elemental tracers. This allows isotopic tracers to be analyzed with classical unmixing models and combined with elemental tracers to potentially increase discriminant capacity.

#### Unmixing and Results

  * **Unmix**: This function assesses the relative contribution of potential sediment sources to each sediment mixture using a mass balance approach. It supports both unconstrained and constrained optimization. The output is a data frame with the relative contributions of sediment sources to each sediment mixture across all iterations.
  * **plot\_results**: This function generates a plot showing the relative contribution of sediment sources to each mixture. It can use either violin charts or density plots.

-----

### Installation

You can install the stable release of `fingerPro` directly from CRAN.

```r
# Install fingerPro from CRAN
install.packages("fingerPro")

# Load the package
library('fingerPro')
```

-----

### Getting Started

Before running any of the package's functions, it's crucial to ensure your data is correctly formatted. The `check_database()` function automatically verifies the integrity of your dataset and infers its type ("raw", "averaged", or "isotopic") based on its column names.

#### Database Formats

The package supports four main database formats, each with specific column requirements:

  * **'raw' format**: Contains individual measurements for scalar tracers. It must have columns for **ID**, **samples**, and **tracer1, tracer2, ...**.
  * **'isotopic raw' format**: Contains individual measurements for isotopic tracers, which require both ratio and content data. It must have columns for **ID**, **samples**, **ratio1, ratio2, ...**, and **cont\_ratio1, cont\_ratio2, ...**.
  * **'averaged' format**: Contains statistical summaries of the scalar tracer data. It must have columns for **ID**, **samples**, **mean\_tracer1, mean\_tracer2, ...**, **sd\_tracer1, sd\_tracer2, ...**, and **n**.
  * **'isotopic averaged' format**: Contains statistical summaries for isotopic tracers. It must have columns for **ID**, **samples**, **mean\_ratio1, mean\_ratio2, ...**, **mean\_cont\_ratio1, mean\_cont\_ratio2, ...**, **sd\_ratio1, sd\_ratio2, ...**, **sd\_cont\_ratio1, sd\_cont\_ratio2, ...**, and **n**.

#### Example Workflow

Here is a typical workflow for using `fingerPro` to analyze a sediment dataset:

1.  **Load and Verify Data**

    ```r
    # Load the example dataset for a 3-source mixing problem.
    data <- read.csv(system.file("extdata", "example_geo_3s_raw.csv", package = "fingerPro"))

    # Verify the structure and integrity of the loaded dataset.
    check_database(data)
    ```

2.  **Exploratory Analysis & Tracer Pre-selection**

      * Visualize tracer distributions with a box plot: `box_plot(data)`.
      * Examine tracer relationships with a correlation plot: `correlation_plot(data)`.
      * Use statistical tests to pre-select discriminant tracers:
          * **Kruskal-Wallis Test**: `KW_test(data)`
          * **Discriminant Function Analysis (DFA)**: `DFA_test(data)`

3.  **Advanced Tracer Selection**

      * Calculate the **Consensus Ranking (CR)** score for each tracer to identify conservative ones: `output_CR <- CR(data, debates = 1000)`.
      * Identify the most discriminant minimal tracer combinations using `CTS_seeds()`: `output_CTS_seeds <- CTS_seeds(data, iter = 1000)`.
      * Perform an **Individual Tracer Analysis (ITA)** to assess tracer behavior in different contexts: `output_ITA <- individual_tracer_analysis(data)`.
      * Based on the ITA results, calculate the **Conservativeness Index (CI)**: `output_CI <- CI(output_ITA)`.

4.  **Final Tracer Selection & Unmixing**

      * Merge the results from the CTS, CR, and CI analyses to create a summary data frame.
      * Filter the summary data to select only the most robust tracers, for instance, by setting criteria for a low CTS error (\< 0.05) and a high CR score (\> 80).
      * Use the `select_tracers()` function to create a new data frame with only the selected tracers: `data <- select_tracers(data, output_data_TracerSelection[, 1])`.
      * Run the unmixing model: `output_UNMIX <- unmix(data)`.

5.  **Visualize and Save Results**

      * Visualize the source contributions using violin or density plots: `plot_results(output_UNMIX)`.
      * Save your results to the workspace with `write_results(output_UNMIX)`.

-----

##  :point_right: Frequently Asked Questions :arrow_right: [Check them!](FaQ.md)

-----

## Contributing and feedback
This software has been improved by the questions, suggestions, and bug reports of the user community. If you have any comments, please use the [Issues](https://github.com/eead-csic-eesa/fingerPro/issues) page or report them to fingerpro@eead.csic.es.

-----

## Citing fingerPro

Latorre, B., Lizaga, I., Gaspar, L., & Navas, A, 2025. Evaluating the Impact of High Source Variability and Extreme Contributing Sources on Sediment Fingerprinting Models. *Water Resources Management*, *1-15*. https://doi.org/10.1007/s11269-025-04169-8

Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2022. Combined use of geochemistry and compound-specific stable isotopes for sediment fingerprinting and tracing. *Science of The Total Environment* 832, 154834. https://doi.org/10.1016/j.scitotenv.2022.154834

Latorre, B., Lizaga, I., Gaspar, L., Navas, A., 2021. A novel method for analysing consistency and unravelling multiple solutions in sediment fingerprinting. *Science of The Total Environment* 789, 147804. https://doi.org/10.1016/j.scitotenv.2021.147804

Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020a. Consensus ranking as a method to identify non-conservative and dissenting tracers in fingerprinting studies. Science of The Total Environment 720, 137537. https://doi.org/10.1016/j.scitotenv.2020.137537

Lizaga, I., Latorre, B., Gaspar, L., Navas, A., 2020. FingerPro: an R package for tracking the provenance of sediment. *Water Resources Management* 272, 111020. https://doi.org/10.1007/s11269-020-02650-0

-----

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

-----

### 🎥 Video Tutorials
[![Alt text](https://img.youtube.com/vi/LcrM_vLOa_I/0.jpg)](https://www.youtube.com/watch?v=LcrM_vLOa_I)

[![Alt text](https://img.youtube.com/vi/7HwGcRSO2O8/0.jpg)](https://www.youtube.com/watch?v=7HwGcRSO2O8&ab_channel=UnmixingScience)
