#' Run DE on counts based on TFIDF
#' @param counts A dgCMatrix with raw or sequencing depth corrected counts
#' @param clusters A character vector of assigned clusters corresponding
#' to columns in the counts matrix
#' @param expression.cutoff Cutoff to determine whether a gene is expressed or not. Default is NULL
#' which uses the median count
#' @importFrom Matrix sparse.model.matrix
#' @importFrom stats as.formula median p.adjust phyper

#' @export
TFIDFMarkers <- function(counts,
                         clusters,
                         expression.cutoff = NULL) {
  counts <- as(object = counts, Class = "dgCMatrix")

  if (is.null(x = expression.cutoff)) {
    expression.cutoff <- median(x = counts@x) - 1
  }
  # Total cells in which gene is expressed
  ncells <- ncol(x = counts)
  cluster.counts <- as.matrix(x = table(clusters))
  cluster.counts.other <- ncells - cluster.counts

  mymodel.matrix <- sparse.model.matrix(object = ~ 0 + clusters)
  colnames(x = mymodel.matrix) <- gsub(pattern = "clusters", replacement = "", x = colnames(x = mymodel.matrix))

  ncells.expressed.percelltype <- (counts > expression.cutoff) %*% mymodel.matrix
  ncells.expressed.percelltype <- as.matrix(x = ncells.expressed.percelltype)
  ncells.expressed.all <- rowSums(x = ncells.expressed.percelltype)

  ncells.expressed.other <- ncells.expressed.all - ncells.expressed.percelltype

  TF <- sweep(x = ncells.expressed.percelltype, MARGIN = 2, STATS = cluster.counts, FUN = "/")
  TF.other <- sweep(x = ncells.expressed.other, MARGIN = 2, STATS = cluster.counts.other, FUN = "/")

  IDF <- log(ncells / ncells.expressed.all)

  TF.IDF <- TF * IDF

  TF <- as.data.frame(x = TF)
  IDF <- as.data.frame(x = IDF)
  TF.IDF <- as.data.frame(x = TF.IDF)

  pvalues <- lapply(colnames(x = TF), function(cluster) {
    phyper(
      q = ncells.expressed.percelltype[, cluster],
      m = ncells.expressed.all,
      n = ncells - ncells.expressed.all,
      k = cluster.counts[cluster, ],
      lower.tail = FALSE
    )
  })


  qvalues <- lapply(pvalues, FUN = p.adjust, method = "BH")
  qvalues <- do.call(cbind, qvalues)
  colnames(qvalues) <- colnames(x = TF)

  pvalues <- do.call(cbind, pvalues)
  colnames(pvalues) <- colnames(x = TF)

  order.tfidf <- lapply(colnames(x = TF), function(e) {
    order(TF.IDF[, e], decreasing = TRUE)
  })
  order.tfidf2 <- cbind(unlist(order.tfidf, use.names = FALSE), rep(seq_len(length.out = ncol(x = TF)), lengths(order.tfidf)))

  de.genes <- rownames(x = TF)[order.tfidf2[, 1]]
  de.clusters <- colnames(x = TF)[order.tfidf2[, 2]]

  de.out <- data.frame(
    gene = de.genes,
    cluster = de.clusters,
    tf.cluster = TF[order.tfidf2],
    tf.other = TF.other[order.tfidf2],
    tf.global = ncells.expressed.all[order.tfidf2[, 1]] / ncells,
    tfidf = TF.IDF[order.tfidf2],
    p_val = pvalues[order.tfidf2],
    p_val_adj = qvalues[order.tfidf2]
  )

  return(de.out)
}

#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows filter
DEOnevsRest <- function(fit, pval_thresh = 0.05) {
  de_list <- list()
  all_celltypes <- grep(pattern = "^celltype", x = colnames(x = fit$Beta), value = TRUE)
  for (celltype in all_celltypes) {
    celltypes_other <- setdiff(all_celltypes, celltype)
    celltypes_other <- paste0("`", celltypes_other, "`")
    celltypes_other <- paste0("(", paste(celltypes_other, collapse = " + "), ")", " / ", length(celltypes_other))

    contrast_other <- paste(celltypes_other, collapse = paste0("/", length(celltypes_other), "+"))
    contrast <- paste0("`", celltype, "`", " - ", contrast_other)
    message(paste0(celltype, ":", contrast))

    de_results <- glmGamPoi::test_de(fit, contrast = contrast) %>%
      arrange(adj_pval) %>%
      filter(pval < pval_thresh)
    de_results$celltype <- gsub(pattern = "celltype", replacement = "", x = celltype)
    de_list[[celltype]] <- de_results
  }
  de_results <- bind_rows(de_list)
  return(de_results)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @export
DE.Pseudobulk <- function(counts, celltypes, covariate_metadata = NULL, pseudocount = 1,
                          k_replicates = 6, seed = 42, pval_thresh = 0.05, offset_libsize = FALSE) {
  # First pseudobulk based on celltype and covariate column so that it is balanced
  # we assign a replicate numbner to each cell such that it has balanced number of covariates
  message("Aggregating counts")
  pseudobulk.out <- PseudobulkCounts(counts,
    celltypes = celltypes, other_covariates_df = covariate_metadata,
    k_replicates = k_replicates, seed = seed
  )
  counts <- pseudobulk.out$counts + pseudocount
  metadata <- pseudobulk.out$metadata
  pct.expressed <- pseudobulk.out$pct.expressed

  formula_str <- "~ 0 + celltype"
  for (name in colnames(covariate_metadata)) {
    formula_str <- paste0(formula_str, " + ", name)
  }
  message("Fitting GLM on aggregated counts")
  if (offset_libsize) {
    offsets <- log(x = Matrix::colSums(x = counts))
    fit <- glmGamPoi::glm_gp(
      data = counts, design = as.formula(formula_str), col_data = metadata,
      offset = offsets, size_factors = FALSE
    )
  } else {
    fit <- glmGamPoi::glm_gp(data = counts, design = as.formula(formula_str), col_data = metadata)
  }

  message("Estimating contrasts")
  de_results <- DEOnevsRest(fit = fit, pval_thresh = pval_thresh)
  return(de_results)
}


#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows
#' @export
GLMDE.SingleCell <- function(counts, celltypes, covariate_metadata = NULL,
                             k_replicates = 6, seed = 42, pval_thresh = 0.05, offset_libsize = TRUE) {
  # First pseudobulk based on celltype and covariate column so that it is balanced
  # we assign a replicate numbner to each cell such that it has balanced number of covariates
  stopifnot(ncol(counts) == length(celltypes))
  metadata.celltype <- data.frame(row.names = colnames(counts), celltype = celltypes)
  if (!is.null(covariate_metadata)) {
    metadata <- cbind(metadata.celltype, covariate_metadata)
  } else {
    metadata <- metadata.celltype
  }
  formula_str <- "~ 0 + celltype"
  if (!is.null(covariate_metadata)) {
    for (name in colnames(covariate_metadata)) {
      formula_str <- paste0(formula_str, " + ", name)
    }
  }
  message("Fitting GLM")
  if (offset_libsize) {
    message("Offsetting libsize")
    offsets <- log(x = Matrix::colSums(x = counts))
    fit <- glmGamPoi::glm_gp(
      data = counts, design = as.formula(formula_str), col_data = metadata, offset = offsets,
      on_disk = FALSE, size_factors = FALSE
    )
  } else {
    fit <- glmGamPoi::glm_gp(data = counts, design = as.formula(formula_str), col_data = metadata, on_disk = FALSE)
  }

  message("Estimating contrasts")
  de_results <- DEOnevsRest(fit = fit, pval_thresh = pval_thresh)
  return(de_results)
}

#' @importFrom magrittr %>%
#' @importFrom dplyr bind_rows filter
#' @export
FilterMarkers <- function(markers, pct.expr, pct.1 = 0.1) {
  ct_list <- list()
  for (ct in unique(markers$celltype)) {
    pct.expr.subset <- pct.expr %>%
      filter(celltype == ct) %>%
      filter(mean_pct > pct.1)
    ct_subset <- markers %>%
      filter(celltype == ct) %>%
      filter(name %in% pct.expr.subset$gene)
    ct_list[[ct]] <- ct_subset
  }
  return(bind_rows(ct_list))
}
