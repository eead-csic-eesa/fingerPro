% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PCAPlot.R
\name{PCAPlot}
\alias{PCAPlot}
\title{Principal component analysis chart}
\usage{
PCAPlot(data, components = c(1:2))
}
\arguments{
\item{data}{Data frame containing source and mixtures data}

\item{components}{Numeric vector containing the index of the two principal components in the chart}
}
\description{
The function performs a principal components analysis on the given data matrix and displays a biplot using vqv.ggbiplot package of the results for each different source to help the user in the decision.
}
