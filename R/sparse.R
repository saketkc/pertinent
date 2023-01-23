#' Sweep for sparse matrices
#' @export
SweepSparse <- function(x,
                        MARGIN,
                        STATS,
                        FUN = "/") {
  if (!inherits(x = x, what = "dgCMatrix")) {
    stop("input needs to be dgCMatrix")
  }
  fun <- match.fun(FUN)
  if (MARGIN == 1) {
    idx <- x@i + 1
    x@x <- fun(x@x, STATS[idx])
  } else if (MARGIN == 2) {
    x <- as(x, "RsparseMatrix")
    idx <- x@j + 1
    x@x <- fun(x@x, STATS[idx])
    x <- as(x, "CsparseMatrix")
  }
  return(x)
}

#' Calculate ranks on sparse matrix
#' @export
SparsifiedRanks <- function(X) {
  if (class(X)[1] != "dgCMatrix") {
    X <- as(object = X, Class = "dgCMatrix")
  }
  non_zeros_per_col <- diff(x = X@p)
  n_zeros_per_col <- nrow(x = X) - non_zeros_per_col
  offsets <- (n_zeros_per_col - 1) / 2
  x <- X@x
  ## split entries to columns
  col_lst <- split(x = x, f = rep.int(1:ncol(X), non_zeros_per_col))
  ## calculate sparsified ranks and do shifting
  sparsified_ranks <- unlist(x = lapply(
    X = seq_along(col_lst),
    FUN = function(i) rank(x = col_lst[[i]]) + offsets[i]
  ))
  ## Create template rank matrix
  X.ranks <- X
  X.ranks@x <- sparsified_ranks
  return(X.ranks)
}


#' Calculate spearman correlation on sparse matrix(ces)
#' @export
SparseSpearmanCor <- function(X, Y = NULL, cov = FALSE) {
  # Get sparsified ranks
  rankX <- SparsifiedRanks(X)
  if (is.null(Y)) {
    # Calculate pearson correlation on rank matrices
    return(corSparse(X = rankX, cov = cov))
  }
  rankY <- SparsifiedRanks(Y)
  return(corSparse(X = rankX, Y = rankY, cov = cov))
}

#' Convert matrix to sparse matrix
#' @export
make.sparse <- function(mat, Class = "CsparseMatrix") {
  # mat <- as(object = mat, Class = "Matrix")
  return(as(object = as(object = as(object = mat, Class = "dMatrix"), Class = "generalMatrix"), Class = Class))
}

#' Convert a sparse matrix to a dataframe long form
#' @export
SparseMatrixToDF <- function(mat) {
  mat <- make.sparse(mat = mat, Class = "TsparseMatrix")
  df <- data.frame(row = rownames(x = mat)[mat@i + 1], column = colnames(x = mat)[mat@j + 1], value = mat@x)
  return(df)
}

#' Convert NAs with zero
#' @export
NAToZero <- function(mat) {
  mat <- replace(mat, is.na(x = mat), 0)
  return(mat)
}


#' Convert mean and variance calculated on non-zero entries into
#' mean and variance if zeroes were to be included. This method is
#' useful when mean and variance are calculated for features in a data
#' frame where the zero entries have been dropped such as when using the
#' \code{SparseMatrixToDF} method.
#' @param nzmean Mean of non-zero entries
#' @param nzvar Mean of non-zero entries
#' @param n Number of non-zero entries
#' @param z Number of zero entries
#' @export
NZstatsToZstats <- function(nzmean, nzvariance, n, z) {
  mean_with_zeros <- n * nzmean / (n + z)
  var_with_zeros <- n * (nzvariance + (nzmean)^2) / (n + z) - mean_with_zeros^2
  return(data.frame(mean = mean_with_zeros, variance = var_with_zeros))
}


#' Merge and add multiple sparse matrices.
#'
#' Assuming the rownames are shared across all matrices, this method creates a collated
#' matrix such that values for shared column names are added while unique columns
#' for each matrices are retained.
#' @param matrix.list A list of sparce matrices
#' @export
AddMergeSparse <- function(matrix.list) {
  mtx1 <- matrix.list[[1]]

  for (i in 2:length(matrix.list)) {
    mtx2 <- matrix.list[[i]]
    mtx1.rownames <- rownames(mtx1)
    mtx2.rownames <- rownames(mtx2)

    stopifnot(mtx1.rownames == mtx2.rownames)

    missing.2 <- colnames(mtx1)[!colnames(mtx1) %in% colnames(mtx2)]
    missing.1 <- colnames(mtx2)[!colnames(mtx2) %in% colnames(mtx1)]

    missing.1.named <- as.vector(x = numeric(length(missing.1)), mode = "list")
    names(missing.1.named) <- missing.1
    missing.2.named <- as.vector(x = numeric(length(missing.2)), mode = "list")
    names(missing.2.named) <- missing.2

    mtx1.expanded <- Reduce(cbind, c(mtx1, missing.1.named))
    mtx2.expanded <- Reduce(cbind, c(mtx2, missing.2.named))
    if (length(missing.1.named) > 0) {
      expanded.length <- ncol(mtx1.expanded) - length(missing.1.named) + 1
      colnames(mtx1.expanded)[expanded.length:ncol(mtx1.expanded)] <- names(missing.1.named)
    }
    if (length(missing.2.named) > 0) {
      expanded.length <- ncol(mtx2.expanded) - length(missing.2.named) + 1
      colnames(mtx2.expanded)[expanded.length:ncol(mtx2.expanded)] <- names(missing.2.named)
    }

    mtx1.expanded <- mtx1.expanded[, sort(colnames(mtx1.expanded))]
    mtx2.expanded <- mtx2.expanded[, sort(colnames(mtx2.expanded))]
    stopifnot(colnames(mtx1.expanded) == colnames(mtx1.expanded))

    mtx1 <- mtx1.expanded + mtx2.expanded
  }
  return(mtx1)
}
