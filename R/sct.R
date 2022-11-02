#' SCT's reconstruction error
#' @param seu Seurat object
#' @export

SCTReconstructionError <- function(seu, dims = 1:30, nfeatures = 3000, idents = c("seurat_clusters"),
                                   outdir = NULL, overwrite = TRUE) {
  sct <- SCTransform(seu, ncells = 5000, variable.features.n = nfeatures, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = dims, verbose = FALSE) %>%
    FindNeighbors(
      dims = dims,
      verbose = FALSE
    ) %>%
    FindClusters(verbose = FALSE)
}



#' @export

fit_phi <- function(mean, variance) {
  # fits a linear model on variance ~ mu + phi*mu^2 to estimate phi (which is
  # 1/theta)
  var_res <- variance - mean
  fit <- lm(var_res ~ 0 + I(mean^2))
  phi <- fit$coefficients[[1]]
  return(list(phi = phi, fit = (mean + phi * mean^2)))
}

#' @import sparseMatrixStats
#' @export

MeanVarFit <- function(counts) {
  means <- sparseMatrixStats::rowMeans2(counts)
  variance <- sparseMatrixStats::rowVars(counts)
  df <- data.frame(mean = means, variance = variance, gene = rownames(counts))
  fit <- fit_phi(means, variance)
  phi <- fit[["phi"]]

  df$fit <- fit[["fit"]]
  df$phi <- phi
  df$residual <- df$variance - df$fit
  return(df)
}




#' @import sparseMatrixStats
#' @export

MeanVarPlot <- function(counts, annotate = FALSE, logxy = F) {
  means <- sparseMatrixStats::rowMeans2(counts)
  variance <- sparseMatrixStats::rowVars(counts)
  df <- data.frame(mean = means, variance = variance, gene = rownames(counts))
  fit <- fit_phi(means, variance)
  phi <- fit[["phi"]]

  df$fit <- fit[["fit"]]

  # this is (mu+mu^2*phi)/mu-1
  genewise_phi <- (variance / means - 1) / means

  compare_phi <- genewise_phi / fit[["phi"]]
  compare_phi[is.na(compare_phi)] <- 0

  index <- genewise_phi > quantile(genewise_phi, 0.999, na.rm = T)[[1]]
  index[is.na(index)] <- FALSE


  df$label <- ""
  df$label[index] <- df$gene[index]

  p <- ggplot(df, aes(x = mean, y = variance, color = "d", label = label)) +
    # geom_point(
    #  alpha = 0.5,
    #  shape = 16, color = "black"
    # )
    scattermore::geom_scattermore(alpha = 0.5)
  if (annotate) {
    p <- p + geom_text_repel()
  }
  p <- p + geom_line(aes(x = mean, y = mean, color = "a"), size = 1.2, show.legend = T) +
    geom_line(aes(x = mean, y = fit, color = "b"), size = 1.2, show.legend = T) +
    xlab("Mean") + ylab("Variance") + ggtitle("Mean-var")
  p <- p + scale_color_manual(name = "", values = c(
    a = "red", b = "blue",
    d = "black"
  ), labels = c("Poisson-fit", paste0("NB-fit (phi=", round(
    phi,
    1
  ), ")"), "Gene")) + theme(legend.position = "bottom")
  if (logxy) {
    p <- p + scale_x_log10() +
      scale_y_log10()
  }
  return(p)
}

#' @import sparseMatrixStats
#' @export

DropoutPlot <- function(counts, annotate = F) {
  n_cells <- dim(counts)[2]
  means <- sparseMatrixStats::rowMeans2(counts)
  variance <- sparseMatrixStats::rowVars(counts)
  dropouts.df <- data.frame(
    mean = means,
    obs_dropout = apply(counts == 0, 1, sum) / n_cells,
    gene = rownames(counts)
  )
  fit <- fit_phi(means, variance)
  phi <- fit[["phi"]]

  dropouts.df$nb_expected_dropout <- (1 / (phi * means + 1))^(1 / phi)

  dropouts.df$pois_expected_dropout <- exp(-dropouts.df$mean)

  expected_over_observed <- log1p(dropouts.df$obs_dropout) - log1p(dropouts.df$pois_expected_dropout)
  index <- expected_over_observed > 0.2
  dropouts.df$label <- ""
  dropouts.df$label[index] <- dropouts.df$gene[index]

  p_counts <- ggplot(dropouts.df, aes(x = mean, y = obs_dropout, color = "d", label = label)) +
    geom_point(alpha = 0.5)
  # scattermore::geom_scattermore(alpha=0.5, )
  if (annotate) {
    p_counts <- p_counts + geom_text_repel()
  }
  p_counts <- p_counts + geom_line(aes(x = mean, y = pois_expected_dropout, color = "a"), size = 1.2, show.legend = T) +
    geom_line(aes(x = mean, y = nb_expected_dropout, color = "b"), size = 1.2, show.legend = T) +
    xlab("mean") + ylab("Dropout probability") + ggtitle("") +
    scale_color_manual(
      name = "",
      values = c("a" = "red", "b" = "blue", "d" = "black"),
      labels = c("Poisson-fit", paste0("NB-fit (phi=", round(phi, 1), ")"), "Gene")
    ) + theme(legend.position = "bottom")
  return(p_counts)
}

#' @import glmGamPoi
#' @export

fit_glmGamPoiOnlyInterceptSingleGene <- function(umi) {
  model_str <- "~ 1"
  fit <- glmGamPoi::glm_gp(
    data = umi, design = as.formula(gsub("y", "", model_str)),
    size_factors = FALSE
  )
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-04)
  colnames(fit$Beta)[1] <- "b0"
  mean <- mean(umi)
  variance <- var(umi)
  overdispersion.factor <- variance / mean
  cv <- sqrt(variance) / mean
  fit.df <- cbind(
    fit$theta, fit$Beta, fit$overdispersions, mean, variance, overdispersion.factor,
    cv
  )
  colnames(fit.df)[1] <- "theta"
  colnames(fit.df)[3] <- "overdispersion"
  colnames(fit.df)[4] <- "mean"
  colnames(fit.df)[5] <- "variance"
  colnames(fit.df)[6] <- "overdispersion.factor"
  colnames(fit.df)[7] <- "cv"
  return(as.data.frame(fit.df))
}

#' @import glmGamPoi
#' @export

getDispersionsOnlyIntercept <- function(counts) {
  gene <- rownames(counts)[1]
  umi <- counts[gene, ]
  fit <- fit_glmGamPoiOnlyInterceptSingleGene(umi)
  fit$gene <- gene
  master.df <- fit
  for (gene in rownames(counts)[2:length(rownames(counts))]) {
    umi <- counts[gene, ]
    fit <- data.frame(fit_glmGamPoiOnlyInterceptSingleGene(umi))
    fit$gene <- gene
    master.df <- rbind(master.df, fit)
  }
  return(master.df)
}

#' @import glmGamPoi
#' @export

fit_glmGamPoiWithLibsizeSingleGene <- function(umi, data) {
  model_str <- "~ log_umi"
  fit <- glmGamPoi::glm_gp(
    data = umi, design = as.formula(gsub("y", "", model_str)),
    col_data = data, size_factors = FALSE
  )
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-04)
  colnames(fit$Beta)[1] <- "b0"
  colnames(fit$Beta)[2] <- "b1"
  fit.df <- cbind(fit$theta, fit$Beta, fit$overdispersions)
  colnames(fit.df)[1] <- "theta"
  colnames(fit.df)[4] <- "overdispersion"
  # means <- mean(umi) variance <- var(umi)

  # fit <- fit_phi(means, variance) phi <- (variance - means)/(means^2)#
  # fit[['phi']] fit.df$phi <- phi
  return(as.data.frame(fit.df))
}


#' @import glmGamPoi
#' @export

getDispersionsWithLibsize <- function(counts) {
  data <- data.frame(log_umi = log10(sparseMatrixStats::colSums2(counts)))
  gene <- rownames(counts)[1]
  umi <- counts[gene, ]
  fit <- fit_glmGamPoiWithLibsizeSingleGene(umi, data)
  cat("\n")
  fit$gene <- gene
  master.df <- fit
  for (gene in rownames(counts)[2:length(rownames(counts))]) {
    umi <- counts[gene, ]
    fit <- data.frame(fit_glmGamPoiWithLibsizeSingleGene(umi, data))
    fit$gene <- gene
    master.df <- rbind(master.df, fit)
  }
  return(master.df)
}
