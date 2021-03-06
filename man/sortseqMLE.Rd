% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE_functions.R
\name{sortseqMLE}
\alias{sortseqMLE}
\title{Maximum likelihood estimation of mu and sigma using sort-seq read count data}
\usage{
sortseqMLE(gates, read.counts, sorted.cells, eps = 0)
}
\arguments{
\item{gates}{List of integers indicating the gate boundaries}

\item{read.counts}{Dataframe with rows corresponding to variants and columns
to gates with elements indicating read counts}

\item{sorted.cells}{List of integers indicating the total number of cells 
sorted into each gate}

\item{eps}{Integer describing the error rate of sorting and sequencing.}
}
\value{
Dataframe with columns indicating the estimated cells for each 
variant sorted into each bin, the total number of cells sorted, and the MLE 
mu and sigma for each variant.
}
\description{
This function is used to obtain maximum likelihood estimates of log-normal 
paramaters mu and sigma for each variant in a sort-seq experiment. The read count
data for each variantis used with the parameters describing sorting gates and
total cells sorted to mininize the log-likelihood function.
}
