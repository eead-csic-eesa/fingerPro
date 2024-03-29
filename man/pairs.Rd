% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairs.R
\name{pairs}
\alias{pairs}
\title{Pairs function}
\usage{
pairs(source, mixture, iter = 1000, seed = 123456)
}
\arguments{
\item{source}{Data frame containing the sediment sources from a dataset}

\item{mixture}{Data frame containing one of the dataset mixtures}

\item{iter}{Number of iteration for each tracer}

\item{seed}{Seed for the random number generator}
}
\value{
Data frame containing all the possible pairs combination from your dataset, their system of equation solution, the consistency and the discriminant capacity.
}
\description{
Extract all the possible selections formed by sets of n ??? 1 tracers and solve them by using standard methods
}
