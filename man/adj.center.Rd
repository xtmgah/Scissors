% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adj.center.R
\name{adj.center}
\alias{adj.center}
\title{Provide a center value adjusted by a lower bound}
\usage{
adj.center(x, average = "mean", trim = 0.1, adjval = NULL)
}
\arguments{
\item{x}{an input vector, data values.}

\item{average}{type of the center being calculated, with choices "median" and "mean".
Default is "median".}

\item{adjval}{the lower bound above which the values will be included to calculate the center.}
}
\description{
This function calculates a center of the values above a given bound.
}
