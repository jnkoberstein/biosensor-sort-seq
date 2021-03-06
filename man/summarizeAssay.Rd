% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE_functions.R
\name{summarizeAssay}
\alias{summarizeAssay}
\title{Combine parameter estimates for both F0 and Fl datasets and calculate the 
dynamic range for each variant}
\usage{
summarizeAssay(F0, Fl)
}
\arguments{
\item{F0}{Dataframe containing maximum likelihood estimates for the 
ligand-free sample.}

\item{Fl}{Dataframe containing maximum likelihood estimates for the 
ligand-bound sample.}
}
\value{
Dataframe combining the two dataframes with an additional column
with the calculated dynamic range.
}
\description{
This function is used to combine the maximum likelihood estimates for the 
two datasets F0 and Fl in order to calculate the dynamic range for each 
variant.
}
