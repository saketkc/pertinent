#' Estimate phi based on a linear fit to mean-variance relationship
#' @importFrom MASS rlm
#' @export
EstimateNBPhi <- function(mean, variance, robust = TRUE) {
  # fits a linear model on variance ~ mu + phi*mu^2 to estimate phi (which is
  # 1/theta)
  var_res <- variance - mean
  if (robust) {
    fit <- suppressWarnings(expr = rlm(formula = var_res ~ 0 + I(mean^2), init = "lts"))
  } else {
    fit <- lm(formula = var_res ~ 0 + I(mean^2))
  }
  phi <- fit$coefficients[[1]]
  return(data.frame(phi = phi, fit = (mean + phi * mean^2)))
}

#' Perform mean-variance calculation on the count matrix.
#' @importFrom sparseMatrixStats rowMeans2 rowVars
#' @export
MeanVarFit <- function(counts = NULL, means = NULL, variance = NULL) {
  if (is.null(x = counts) & ((is.null(x = means) | is.null(x = variance)))) {
    stop("Need to specify one of counts or (means, variance).")
  }
  if (!is.null(x = counts)) {
    means <- rowMeans2(x = counts)
    variance <- rowVars(x = counts)
  }

  df <- data.frame(mean = means, variance = variance, gene = rownames(counts))
  fit <- EstimateNBPhi(means, variance)
  df$fit <- fit$fit
  df$phi <- fit$phi
  df$residual <- df$variance - df$fit
  return(df)
}


#' Plot mean variance relationship for counts or dataframe
#' @importFrom sparseMatrixStats rowMeans2 rowVars
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle xlab ylab theme scale_color_manual
#' @export
MeanVarPlot <- function(counts = NULL, df = NULL, annotate = FALSE, logxy = T) {
  if (is.null(x = df)) {
    means <- rowMeans2(x = counts)
    variance <- rowVars(x = counts)
    df <- data.frame(mean = means, variance = variance, gene = rownames(counts))
  }
  fit <- EstimateNBPhi(df$mean, df$variance)
  phi <- fit$phi
  df$fit <- fit$fit

  # this is (mu+mu^2*phi)/mu-1
  genewise_phi <- (variance / means - 1) / means
  compare_phi <- genewise_phi / fit$phi
  compare_phi[is.na(compare_phi)] <- 0
  index <- genewise_phi > quantile(genewise_phi, 0.999, na.rm = T)[[1]]
  index[is.na(index)] <- FALSE

  df$label <- ""
  df$label[index] <- df$gene[index]

  p <- ggplot(df, aes(x = mean, y = variance, color = "d", label = label)) +
    geom_point(
      alpha = 0.5,
      shape = 16,
      color = "black"
    )

  if (annotate) {
    p <- p + geom_text_repel()
  }
  p <- p + geom_line(aes(x = mean, y = mean, color = "a"), size = 1.2, show.legend = T) +
    geom_line(aes(x = mean, y = fit, color = "b"), size = 1.2, show.legend = T) +
    xlab("Mean") +
    ylab("Variance") +
    ggtitle("")
  p <- p + scale_color_manual(
    name = "", values = c(a = "#1A85FF", b = "#D41159", d = "black"),
    labels = c("Poisson-fit", paste0("NB-fit (phi=", round(phi, 1), ")"), "Gene")
  ) +
    theme(legend.position = "bottom")
  if (logxy) {
    p <- p + scale_x_log10() +
      scale_y_log10()
  }
  return(p)
}

#' Plot expected vs observed dropouts (zero observation)
#' @importFrom sparseMatrixStats rowMeans2 rowVars
#' @importFrom ggplot2 ggplot aes geom_point geom_line ggtitle xlab ylab theme scale_color_manual
#' @export
DropoutPlot <- function(counts = NULL, annotate = FALSE) {
  n_cells <- dim(counts)[2]
  means <- rowMeans2(x = counts)
  variance <- rowVars(x = counts)
  dropouts.df <- data.frame(
    mean = means,
    obs_dropout = apply(counts == 0, 1, sum) / n_cells,
    gene = rownames(counts)
  )
  fit <- EstimateNBPhi(means, variance)
  phi <- fit$phi

  dropouts.df$nb_expected_dropout <- (1 / (phi * means + 1))^(1 / phi)
  dropouts.df$pois_expected_dropout <- exp(-dropouts.df$mean)

  expected_over_observed <- log1p(x = dropouts.df$obs_dropout) - log1p(x = dropouts.df$pois_expected_dropout)
  index <- expected_over_observed > 0.2
  dropouts.df$label <- ""
  dropouts.df$label[index] <- dropouts.df$gene[index]

  p_counts <- ggplot(dropouts.df, aes(x = mean, y = obs_dropout, color = "d", label = label)) +
    geom_point(
      alpha = 0.5,
      shape = 16,
      color = "black"
    )
  if (annotate) {
    p_counts <- p_counts + geom_text_repel()
  }
  p_counts <- p_counts +
    geom_line(aes(x = mean, y = pois_expected_dropout, color = "a"), size = 1.2, show.legend = T) +
    geom_line(aes(x = mean, y = nb_expected_dropout, color = "b"), size = 1.2, show.legend = T) +
    xlab("mean") +
    ylab("Dropout probability") +
    ggtitle("") +
    scale_color_manual(
      name = "",
      values = c(a = "#1A85FF", b = "#D41159", "d" = "black"),
      labels = c("Poisson-fit", paste0("NB-fit (phi=", round(phi, 1), ")"), "Gene")
    ) + theme(legend.position = "bottom")
  return(p_counts)
}

#' @import glmGamPoi
#' @export

fit_glmGamPoiOnlyInterceptSingleGene <- function(umi) {
  model_str <- "~ 1"
  fit <- glm_gp(
    data = umi, design = as.formula(gsub("y", "", model_str)),
    size_factors = FALSE
  )
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-04)
  colnames(x = fit$Beta)[1] <- "b0"
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
  fit <- glm_gp(
    data = umi, design = as.formula(gsub("y", "", model_str)),
    col_data = data, size_factors = FALSE
  )
  fit$theta <- pmin(1 / fit$overdispersions, rowMeans(fit$Mu) / 1e-04)
  colnames(fit$Beta)[1] <- "b0"
  colnames(fit$Beta)[2] <- "b1"
  fit.df <- cbind(fit$theta, fit$Beta, fit$overdispersions)
  colnames(fit.df)[1] <- "theta"
  colnames(fit.df)[4] <- "overdispersion"

  return(as.data.frame(x = fit.df))
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
