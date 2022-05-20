################################################################################
# Code example for the YouTube video Unmixing science 2. The FingerPro package
################################################################################

# GitHub repository: https://github.com/eead-csic-eesa/fingerPro
# Citing FingerPro and its tools: https://github.com/eead-csic-eesa/fingerPro#citing-fingerpro-and-its-tools
# ResearchGate: https://www.researchgate.net/profile/Ivan-Lizaga

# The creation of this code was supported by the Research Foundation - Flanders, mandate 12V8622N at 
# Ghent University/Department of Green Chemistry and Technology/ISOFYS group
################################################################################

# To ensure you're working with the last version of the package install it either from GitHub 
# or download the .zip file from GitHub repository and install it in your computer

# From CRAN (version 1.1 || 20/05/2022)
install.packages("fingerPro")

# Download the fingerPro_1.3.zip file from GitHub on you computer
setwd("C:/your/file/directory")
install.packages('fingerPro_1.3.zip', repos = NULL)

# From GitHub (version 1.3)
devtools::install_github("eead-csic-eesa/fingerPro", ref = "master", force = T)


library(fingerPro)
?fingerPro

# Example stored in the help file
################################################################################
# Created by Ivan Lizaga 28/01/2022

# If you want to use your own data,
# I recommend following the format of the datasets included in the package
# setwd("the directory that contains your dataset")
# data <- read.table('your dataset.csv', header = T, sep = ',')

# Example of the data included in the fingerPro package
# Load the dataset called "catchment" 
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
################################################################################
rm(list = ls())
dev.off()

################################################################################
# Compute the CI method and plot the triangles
################################################################################
# The following database was extracted from Raigani et al. (2019): https://doi.org/10.1016/j.ejrh.2019.100613
sources.file <- system.file("extdata", "Raigani.csv", package = "fingerPro")
data <- read.csv(sources.file)
rm(sources.file)

results_CI <- CI_Method(data, points = 2000, Means = T)

Ternary_diagram(results_CI, tracers = c(1:8), n_row = 1, n_col = 8, Means = T)

################################################################################
# Compute the CR method
################################################################################
mixture <- tail(data, n = 1)
var <- grep("^D", colnames(data))
mixture <- mixture[-c(var)]
row.names(mixture) <- NULL
source <- head(data,-1)
sgeo <- source[, -1]
mgeo <- mixture[, -1]

crgeo <- cr_3s(source=sgeo, mixture=mgeo, maxiter = 1000, seed=1234567)# 2000 it's more recommended
print(crgeo)


################################################################################
# Compute the CTS method
################################################################################
library(dplyr)
# compute pairs
# Be aware that for four sources, instead of pairs of tracers you'll have triplets
pgeo <- pairs(sgeo, mgeo, iter = 1000, seed=1234567) # 2000 it's more recommended
head(pgeo)

#Explore those pairs
sol <- pgeo[pgeo$id=="P Sr",]
ctsgeo <- cts_3s(source=sgeo, mixture=mgeo, sol=c(sol$w1, sol$w2, sol$w3))
ctsgeo <- ctsgeo %>% right_join(crgeo, by=c("tracer"))
ctsgeo <- ctsgeo[ctsgeo$err<0.025 & ctsgeo$score>80,]
head(ctsgeo)

sol <- pgeo[pgeo$id=="Ba Sr",]
ctsgeo <- cts_3s(source=sgeo, mixture=mgeo, sol=c(sol$w1, sol$w2, sol$w3))
ctsgeo <- ctsgeo %>% right_join(crgeo, by=c("tracer"))
ctsgeo <- ctsgeo[ctsgeo$err<0.025 & ctsgeo$score>80,]
head(ctsgeo)


################################################################################
# Unmix using the Consistent solutions
################################################################################
# First we select the tracers suggested by the CTS and CR together with their SD (D*) for the two first pairs
data1 <- select(data,"id","sources","Li", "P", "Sr", "Ca", "DLi", "DP", "DSr", "DCa", "n")
# result_FP_1 <- unmix_R(data1, iter = "medium", Means = T) # New version, coming soon
result_FP_1 <- unmix(data1, samples = 200 ,iter = 300, Means = T)

data1 <- select(data,"id","sources","Ba", "Li", "K", "Sr", "DBa", "DLi", "DK", "DSr", "n")
# result_FP_2 <- unmix_R(data1, iter = "medium", Means = T) # New version, coming soon
result_FP_2 <- unmix(data1, samples = 200 ,iter = 300, Means = T)

library(gridExtra)
grid.arrange(plotResults(result_FP_1, y_high = 1), plotResults(result_FP_2, y_high = 1),
             plotResults(result_FP_1, violin = T), plotResults(result_FP_2, violin = T), ncol=2, nrow =2)

# In the previous example we saw two conservative, consensual and consistent solutions that despite being different had
# similar results. This is not always the case as can be seen in the CTS paper. Thus, you can find totally different
# results for your mixture, mostly depending on the nature of your tracers, your sources, transport processes..... 
# and now that you know of their existence you have the opportunity to discuss about why!!!
# How to choose between them?
# 1) Not everything is black or white and there is not unique or define method for that, 
#     display them in your research and highlight this characteristic, let's move towards a less black box unmixing. 
# 2) Use all previous tools combined with the CI, CR and the triangles together with your expert knowledge, 
#     the unmixing is only the icing on the cake after you understand your data.
# 3) New methods and improvements are being constantly uploaded to the repository, so in case of doubts do not hesitate
#     to contact the developers || lizaga.ivan10@gmail.com