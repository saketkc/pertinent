% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparse.R
\name{AddMergeSparse}
\alias{AddMergeSparse}
\title{Merge and add multiple sparse matrices.}
\usage{
AddMergeSparse(matrix.list)
}
\arguments{
\item{matrix.list}{A list of sparce matrices}
}
\description{
Assuming the rownames are shared across all matrices, this method creates a collated
matrix such that values for shared column names are added while unique columns
for each matrices are retained.
}
