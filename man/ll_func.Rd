% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE_functions.R
\name{ll_func}
\alias{ll_func}
\title{Log-likelihood function for binned sort-seq}
\usage{
ll_func(par, ri, lj, uj, dj, eps)
}
\arguments{
\item{par}{parameters to be optimized}

\item{ri}{read counts for variant i}

\item{lj}{lower boundaries for j gates}

\item{uj}{upper boundaries for j gates}

\item{dj}{proportionality factor describing the total number of cell sorted 
to the total number of read counts in each bin.}

\item{eps}{error rate for sorting}
}
\value{
integer representing function evaluated at tested parameters
}
\description{
Log-likelihood function for binned sort-seq
}
\keyword{internal}
