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
  return(list(tau = tau, gap = gap))
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

#' Project a given matrix along a vector
#' The vector is first standardized
#' @export
ProjectDataLongVector <- function(mtx.genebycell, f.vector) {
  f.vector <- f.vector - mean(x = f.vector)
  f.vector <- f.vector / sqrt(x = sum(f.vector^2))
  projection <- (mtx.genebycell - colSums(x = mtx.genebycell)) %*% f.vector
  return(projection)
}

#' Calculate reconstruction error for a new dataset given feature loadings
#' @importFrom Matrix t
#' @export
PredictPCA <- function(feature.loadings, newdata) {
  # rotation is cell embedding
  # x is gene embedding if using prcomp results
  new.embeddings <- t(newdata) %*% feature.loadings
  colnames(new.embeddings) <- colnames(feature.loadings)
  reconstructed.data <- feature.loadings %*% t(new.embeddings)
  reconstruction.error <- sum((newdata - reconstructed.data)^2) / ncol(newdata)
  return(reconstruction.error)
}

#' Perform PCA and return loadings along with reconstruction error
#' @importFrom irlba irlba prcomp_irlba
#' @importFrom Matrix t
#' @export
#'
ReconErrorPCA <- function(original.data, npcs = 50, weight.by.var = TRUE) {
  pca.results <- irlba(A = t(x = original.data), nv = npcs)
  feature.loadings <- pca.results$v
  sdev <- pca.results$d / sqrt(max(1, ncol(original.data) - 1))
  variance.explained <-
    weight.by.var <- TRUE
  if (weight.by.var) {
    cell.embeddings <- pca.results$u %*% diag(pca.results$d)
  } else {
    cell.embeddings <- pca.results$u
  }

  reconstructed.data <- feature.loadings %*% t(cell.embeddings)
  reconstruction.error <- sum((original.data - reconstructed.data)^2) / ncol(original.data)
  return(list(
    error = reconstruction.error,
    feature.loadings = feature.loadings,
    cell.embeddings = cell.embeddings
  ))
}

#' Perform PCA
#' @importFrom irlba irlba prcomp_irlba
#' @export
DoPCA <- function(mtx.genebycell, npcs = 50,
                  weight.by.var = TRUE,
                  do.scale = FALSE, do.center = TRUE,
                  method = "prcomp_irlba") {
  mtx.genebycell <- Seurat::ScaleData(mtx.genebycell,
    do.scale = do.scale,
    do.center = do.center
  )
  npcs <- min(npcs, min(nrow(mtx.genebycell), ncol(mtx.genebycell)))
  if (method == "irlba") {
    npcs <- min(npcs, min(nrow(mtx.genebycell) - 1, ncol(mtx.genebycell) - 1))
    pca.results <- irlba(A = t(x = mtx.genebycell), nv = npcs)
    feature.loadings <- pca.results$v
    sdev <- pca.results$d / sqrt(max(1, ncol(mtx.genebycell) - 1))
    variance.explained <- sdev^2
    total.variance <- NA
    if (weight.by.var) {
      cell.embeddings <- pca.results$u %*% diag(pca.results$d)
    } else {
      cell.embeddings <- pca.results$u
    }
  } else if (method == "prcomp_irlba") {
    npcs <- min(npcs, min(nrow(mtx.genebycell) - 1, ncol(mtx.genebycell) - 1))
    pca.results <- prcomp_irlba(
      x = t(x = mtx.genebycell), n = npcs,
      center = do.center, scale. = do.scale
    )
    feature.loadings <- pca.results$rotation
    sdev <- pca.results$sdev
    variance.explained <- sdev^2
    total.variance <- pca.results$totalvar
    if (weight.by.var) {
      cell.embeddings <- pca.results$x
    } else {
      cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = mtx.genebycell) - 1))
    }
  } else {
    pca.results <- prcomp(x = t(mtx.genebycell), rank. = npcs, ...)
    feature.loadings <- pca.results$rotation
    sdev <- pca.results$sdev
    variance.explained <- sdev^2
    total.variance <- sum(variance.explained)
    if (weight.by.var) {
      cell.embeddings <- pca.results$x
    } else {
      cell.embeddings <- pca.results$x / (pca.results$sdev[1:npcs] * sqrt(x = ncol(x = mtx.genebycell) - 1))
    }
  }
  return(list(
    variance.explained = variance.explained,
    total.variance = total.variance,
    feature.loadings = feature.loadings,
    cell.embeddings = cell.embeddings
  ))
}


#' Regress out covariates (or batches)
#' @importFrom stats lm
#' @export
RegressOutLM <- function(data, cols = NULL, vars.to.regress) {
  # not optimized for sparse
  data <- as.data.frame(data)
  stopifnot(vars.to.regress %in% colnames(data))
  if (is.null(x = cols)) {
    cols <- colnames(data)
  }
  data.subset <- data[, setdiff(cols, vars.to.regress)]
  covariates <- data[, vars.to.regress, drop = FALSE]
  data.regress <- apply(data.subset, 2, function(y) {
    df <- cbind(y, covariates)
    fit <- lm(y ~ ., data = df)
    fit.residuals <- fit[["residuals"]] + fit[["coefficients"]][1]
    return(as.vector(fit.residuals))
  })
  return(data.regress)
}

#' Variance explained in PCA
#' @export
VarianceExplainedPCA <- function(mtx.genebycell, col_data, factors.of.interest, ...) {
  data.pca <- DoPCA(mtx.genebycell = mtx.genebycell, ...)
  frac.variance.explained <- rep(NA, length(factors.of.interest))
  names(frac.variance.explained) <- factors.of.interest
  for (foi in factors.of.interest) {
    foi.vector <- col_data[, foi]
    if (!is.numeric(foi.vector)) {
      foi.vector <- factor(foi.vector)
      foi.vector <- as.matrix(model.matrix(~ 0 + foi.vector))
      # pca.results <- prcomp(foi.matrix, center = TRUE, scale. = FALSE)
    }
    projected <- ProjectDataLongVector(
      mtx.genebycell,
      foi.vector
    )
    frac.variance.explained[[foi]] <- var(projected) / data.pca$total.variance
  }
  return(frac.variance.explained)
}

