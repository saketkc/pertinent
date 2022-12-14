% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudobulk.R
\name{CreateKPseudoBulk}
\alias{CreateKPseudoBulk}
\title{Create a pseudobulk of input matrix given cell.group}
\usage{
CreateKPseudoBulk(
  counts,
  cell.group,
  k.replicates = 6,
  seed = 42,
  return.seurat = FALSE
)
}
\arguments{
\item{counts}{A sparse matrix}

\item{cell.group}{A vector of assigned group/group (one per cell)}

\item{k.replicates}{Number of replicates}

\item{seed}{Integer seed}

\item{return.seurt}{bool Whether seurat object should be returned}
}
\description{
Create a pseudobulk of input matrix given cell.group
}
\examples{
data("pbmc_small")
counts <- GetAssayData(object = pbmc_small, assay = "RNA", slot = "counts")
cell.group <- pbmc_small$groups
k.replicates <- 3
counts.agg <- CreateKPseudoBulk(counts = counts, cell.group = cell.group, k.replicates = k.replicates)
}
