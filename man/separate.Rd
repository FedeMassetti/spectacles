% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spectra-methods.R
\docType{methods}
\name{separate}
\alias{separate}
\alias{separate.Spectra}
\alias{separate,Spectra-method}
\title{Separates a \code{Spectra*} object into a calibration and a validation set.}
\usage{
\S4method{separate}{Spectra}(obj, calibration)
}
\arguments{
\item{obj}{an object inheriting from class \code{SpectraDataFrame}}

\item{calibration}{The fraction of the dataset to be put in the calibration
set}

\item{...}{Ignored}
}
\value{
An list with two \code{SpectraDataFrame} objects, one for the
calibration, and the other for the validation.
}
\description{
Separates a \code{Spectra*} object into a calibration and a validation set.
}
\examples{

# Loading example data
data(australia)
spectra(australia) <- sr_no ~ ... ~ 350:2500

l <- separate(australia, calibration=0.7)
# The result is a list of two Spectra* objects
str(l)
lapply(l, nrow)
}
\author{
Pierre Roudier \url{pierre.roudier@gmail.com}
}
