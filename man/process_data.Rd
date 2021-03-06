% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_data.R
\name{process_data}
\alias{process_data}
\title{Processing data}
\usage{
process_data(rawdata, exon, input.type = "whole_intron",
  output.type = "part_intron", intron.len = NULL,
  logshift.val = NULL, param.grid = NULL, average = "mean",
  trim = 0.1, adjval = NULL, center.type = 1, center.adjval = NULL,
  loop = TRUE, smoothness = 0.7, draw.plot = FALSE,
  plot.main = "Gene", weigh = FALSE, weight.method = "mean",
  max.weight = 2)
}
\description{
Processing data
}
