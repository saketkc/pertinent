#' Calculate Tau index to quantify the specificity of a gene to a cluster
#' @param x A vector of gene's normalized expression
#' @returns A scalar between 0-1 where 1 indicates ubiquitous expression.
#' @export
TauIndex <- function(x) {
  # tau = \sum_{i=1}^n (1-xhat_i)/(n-1) where xhat = x_i/max(x_i)
  # See: https://academic.oup.com/bib/article/18/2/205/2562739
  xhat <- x / max(x, na.rm = T)
  # The ‘gap’ index for x is the maximum difference between the two neighboring values in the sorted vector.
  # If the same ‘gap’ is found more than once in a profile, the first gap, between the smaller neighboring
  # values with that is considered taken. The ‘gap’ can be used to convert expression profiles into binary form.
  index <- order(x)
  xsorted <- x[index]
  gaps <- c(0, diff(x = xsorted))
  gap <- xsorted[which(gaps == min(gaps))]
  n <- length(x)
  tau <- sum(1 - xhat) / (n - 1)
  return(list(tau=tau, gap=gap))
}

#' Calculate gini index to quantify the specificity of a gene to a cluster
#' @param x A vector of gene's normalized expression
#' @param weights A vector of weights for observations (default is 1)
#' @returns A scalar between 0-1 where 1 indicates ubiquitous expression.
#' @export
GiniIndex <- function(x, weights = NULL) {
  # See: https://core.ac.uk/download/pdf/41339501.pdf

  if (is.null(x = weights)) {
    weights <- rep(x = 1, length(x = x))
  }
  weights <- weights / sum(weights)
  index <- order(x)
  x <- x[index]
  weights <- weights[index]

  Fhat <- weights / 2 + c(0, head(x = cumsum(weights), -1))

  xbar <- sum(weights * x)
  Fbar <- sum(weights * Fhat)

  gini <- 2 / xbar * sum(weights * (x - xbar) * (Fhat - Fbar))
  return(gini)
}
