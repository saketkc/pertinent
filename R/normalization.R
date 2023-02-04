#' Proportional filtering normalization
#' @importFrom Matrix colSums
#' @export
NormalizePF <- function(counts) {
  total_umi <- colSums(counts)
  pf <- SweepSparse(counts, MARGIN = 2, STATS = total_umi / mean(total_umi), FUN = "/")
  return(pf)
}

#' Log1p Proportional filtering normalization
#' @importFrom Matrix colSums
#' @export
NormalizeLog1pPF <- function(counts) {
  total_umi <- colSums(counts)
  pf <- SweepSparse(counts, MARGIN = 2, STATS = total_umi / mean(total_umi), FUN = "/")
  log1p_pf <- log1p(pf)
  return(log1p_pf)
}

#' PFLog1pPF normalization
#' @importFrom Matrix colSums
#' @export
NormalizePFLog1pPF <- function(counts) {
  total_umi <- colSums(counts)
  pf <- SweepSparse(counts, MARGIN = 2, STATS = total_umi / mean(total_umi), FUN = "/")
  log1p_pf <- log1p(pf)

  pf_log1p_pf <- NormalizePF(log1p_pf)
  return(pf_log1p_pf)
}

#' Sqrt normalization
#' @importFrom Matrix colSums
#' @export
NormalizeSqrt <- function(counts) {
  return(sqrt(counts))
}

#' Log1p normalization
#' @importFrom Matrix colSums
#' @export
NormalizeLog1p <- function(counts) {
  return(log1p(counts))
}

#' Log1p CP10k normalization
#' @importFrom Matrix colSums
#' @export
NormalizeLog1pCP10k <- function(counts) {
  total_umi <- colSums(counts)
  target_sum <- 10000
  counts.norm <- SweepSparse(counts, MARGIN = 2, STATS = total_umi / target_sum, FUN = "/")
  return(log1p(counts.norm))
}

#' Log1p CPM
#' @importFrom Matrix colSums
#' @export
NormalizeLog1pCPM <- function(counts) {
  total_umi <- colSums(counts)
  target_sum <- 1000000
  counts.norm <- SweepSparse(counts, MARGIN = 2, STATS = total_umi / target_sum, FUN = "/")
  return(log1p(counts.norm))
}

#' Reverse log normalization
#' @importFrom Matrix colSums
#' @export
ReverseLogNorm <- function(lognorm_data) {
  normdata <- as(object = expm1(lognorm_data), Class = "dgCMatrix")
  counts.pseudo <- SweepSparse(normdata,
    MARGIN = 2,
    STATS = apply(normdata, FUN = function(col) {
      min(col[col > 0])
    }, MARGIN = 2),
    FUN = "/"
  )
  return(counts.pseudo)
}
