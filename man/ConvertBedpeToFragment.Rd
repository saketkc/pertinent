% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atac.R
\name{ConvertBedpeToFragment}
\alias{ConvertBedpeToFragment}
\title{Convert a given input bed pe file to a fragment file
to make it comptabile with downstream scATAC-seq workflows}
\usage{
ConvertBedpeToFragment(
  bedpe.path,
  fragment.path,
  chroms.to.keep = c(paste0("chr", seq(1, 22)), paste0("chr", c("X", "Y")))
)
}
\arguments{
\item{bedpe.path}{Path to input bedpe file (can be gzipped)}

\item{fragment.path}{Path to output fragment file (gz compressed; sorted and tabix indexed)}
}
\description{
Convert a given input bed pe file to a fragment file
to make it comptabile with downstream scATAC-seq workflows
}
