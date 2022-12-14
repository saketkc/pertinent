% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/peaks.R
\name{MergePeaksIterative}
\alias{MergePeaksIterative}
\title{Iteratively merge peaks extended around summits
Originally proposed in Corces and Granja 2018 (10.1126/science.aav1898)}
\usage{
MergePeaksIterative(
  summits.path.list = NULL,
  summits.gr.list = NULL,
  bsgenome = NULL,
  blacklist = NULL,
  extend = 250,
  min.spm = 1,
  ...
)
}
\description{
Iteratively merge peaks extended around summits
Originally proposed in Corces and Granja 2018 (10.1126/science.aav1898)
}
