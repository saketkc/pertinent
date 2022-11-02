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
