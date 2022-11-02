#' Summarize difference in prediction obtained from two prediction methods
#' @param object Seurat object
#' @param predcol.1 Metadata column name for first method's predictions
#' @param predcol.2 Metadata column name for second methods' predictions
#' @param min.diff.cells Return only comparisons where the difference between two methods exceeds this. Default is 3
#' @param diff.only Return statistics only for cases where the predictions differ. Default is false which will return all comparisons
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange desc
#' @export
SummarisePredDiffs <- function(object,
                               predcol.1,
                               predcol.2,
                               min.diff.cells = 3,
                               diff.only = FALSE) {
  preds <- as.data.frame(table(
    object@meta.data[, predcol.1],
    object@meta.data[, predcol.2]
  ))

  preds <- preds %>%
    filter(Freq >= min.diff.cells) %>%
    arrange(Var1, Var2, desc(Freq))
  preds$Var1 <- as.character(preds$Var1)
  preds$Var2 <- as.character(preds$Var2)
  same_celltypes <- preds[preds$Var1 == preds$Var2, ]

  colnames(preds)[1] <- predcol.1
  colnames(preds)[2] <- predcol.2

  if (diff.only) {
    preds <- preds[preds[, predcol.2] != preds[, predcol.1], ] %>% arrange(desc(Freq))
    preds <- preds[preds[, predcol.1] %in% same_celltypes$Var1, ]
    preds <- preds[preds[, predcol.2] %in% same_celltypes$Var1, ]
  }

  return(preds)
}

#' Compare accuracy of prediction of two different annotation methods
#'
#' This method uses consensus prediction between the methods as a ground truth to find markers for
#' each of the two celltypes where the two methods have a consensus and then uses these markers
#' to predict the most likely celltype of cells where the predictions differ
#' @param object Seurat object
#' @param celltype.1 Celltype 1
#' @param celltype.2 Celltype 2
#' @param predcol.1 Metadata column name for first method's predictions
#' @param predcol.2 Metadata column name for second methods' predictions
#' @param predcol.1 Metadata column name for first method's prediction scores
#' @param predcol.2 Metadata column name for second methods' prediction scores
#' @param use_scale_data_for_glm Whether to use scaled data for running logistic regression
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange desc top_n
#' @importFrom Matrix t
#' @importFrom stats glm predict
#' @export
ComparePredictions <- function(object,
                               celltype.1,
                               celltype.2,
                               predcol.1,
                               predcol.2,
                               scorecol.1,
                               scorecol.2,
                               assay = "SCT",
                               slot = "data",
                               use.presto = FALSE,
                               use.fastglm = FALSE,
                               do.plot = FALSE) {
  celltypes_to_use <- c(celltype.1, celltype.2)

  object$concat_pred <- paste0(object@meta.data[, predcol.1], "_", object@meta.data[, predcol.2])

  same_celltype <- c(
    paste0(celltype.1, "_", celltype.1),
    paste0(celltype.2, "_", celltype.2)
  )
  diff_celltype <- c(
    paste0(celltype.1, "_", celltype.2),
    paste0(celltype.2, "_", celltype.1)
  )

  obj.subset <- subset(x = object, concat_pred %in% c(same_celltype, diff_celltype))
  obj.subset$concat_pred <- factor(obj.subset$concat_pred, levels = c(
    same_celltype[1], diff_celltype[1],
    same_celltype[2], diff_celltype[2]
  ))
  Idents(object = obj.subset) <- "concat_pred"


  if (use.presto) {
    DefaultAssay(obj.subset) <- assay
    markers <- SeuratWrappers::RunPresto(obj.subset,
      ident.1 = same_celltype[1], ident.2 = same_celltype[2], assay = assay, slot = slot,
      group.by = "concat_pred", min.pct = 0.1, min.cells.group = 1, verbose = FALSE
    )
  } else {
    markers <- FindMarkers(obj.subset,
      ident.1 = same_celltype[1], ident.2 = same_celltype[2], assay = assay, slot = slot,
      group.by = "concat_pred", min.pct = 0.1, min.cells.group = 1, verbose = FALSE
    )
  }
  markers$gene <- rownames(markers)

  markers <- markers[markers$p_val < 0.05, ]
  features.to.plot <- markers$gene

  obj.same_pred <- subset(obj.subset, concat_pred %in% same_celltype)
  obj.nonsame_pred <- subset(obj.subset, concat_pred %in% diff_celltype)

  # Use only 500 max genes for doing logistic regression
  # markers_shortlist <- markers %>% top_n(n = 500, wt = abs(avg_log2FC))
  markers_pos <- markers %>%
    filter(avg_log2FC > 0) %>%
    top_n(n = 100, wt = abs(avg_log2FC))

  markers_neg <- markers %>%
    filter(avg_log2FC < 0) %>%
    top_n(n = 100, wt = abs(avg_log2FC))

  # markers_shortlist <- markers %>% top_n(n = 500, wt = abs(avg_log2FC))
  markers_shortlist <- rbind(markers_pos, markers_neg)
  genes <- unique(as.character(markers_shortlist$gene))


  #
  train.data <- as.data.frame(x = Matrix::t(x = GetAssayData(object = obj.same_pred, assay = assay, slot = slot)[genes, ]))
  # if (use_scale_data_for_glm) {
  #  train.data <- as.data.frame(Matrix::t(obj.same_pred@assays$RNA@scale.data[genes, ]))
  # } else {
  #  train.data <- as.data.frame(Matrix::t(obj.same_pred@assays$SCT@data[genes, ]))
  # }

  # both columns are same
  train.data$celltype <- factor(obj.same_pred@meta.data[, predcol.2], levels = celltypes_to_use)

  test.data <- as.data.frame(x = Matrix::t(x = GetAssayData(object = obj.nonsame_pred, assay = assay, slot = slot)[genes, ]))

  # if (use_scale_data_for_glm) {
  #  test.data <- as.data.frame(Matrix::t(obj.nonsame_pred@assays$RNA@scale.data[genes, ]))
  # } else {
  #  test.data <- as.data.frame(Matrix::t(obj.nonsame_pred@assays$SCT@data[genes, ]))
  # }

  if (use.fastglm) {
    x <- train.data %>%
      select(-celltype) %>%
      as.matrix()
    y <- as.numeric(train.data$celltype) - 1
    model <- fastglm::fastglm(x, y, family = binomial(link = "logit"))
    preds <- predict(model, newdata = as.matrix(x = test.data), type = "response")
  } else {
    model <- glm(as.formula("celltype ~ ."), family = binomial(link = "logit"), data = train.data)
    preds <- predict(model, newdata = test.data, type = "response")
  }


  test.data$pred.1 <- factor(obj.nonsame_pred@meta.data[, predcol.1], levels = celltypes_to_use)
  test.data$pred.2 <- factor(obj.nonsame_pred@meta.data[, predcol.2], levels = celltypes_to_use)

  test.data$prediction.logistic <- celltypes_to_use[1] # this is level 1
  test.data$prediction.logistic[preds > 0.5] <- celltypes_to_use[2] # this is level 2
  test.data$prediction.logistic <- factor(test.data$prediction.logistic, levels = celltypes_to_use)
  test.data$predscore.logistic <- as.numeric(preds)
  test.data$predscore.1 <- obj.nonsame_pred@meta.data[, scorecol.1]
  test.data$predscore.2 <- obj.nonsame_pred@meta.data[, scorecol.2]

  test.data$cell_barcode <- colnames(obj.nonsame_pred@assays$SCT@data)

  matching.1 <- length(which(test.data$pred.1 == test.data$prediction.logistic))
  matching.2 <- length(which(test.data$pred.2 == test.data$prediction.logistic))

  accuracy.1 <- matching.1 / nrow(test.data)
  accuracy.2 <- matching.2 / nrow(test.data)




  if (do.plot) {
    DefaultAssay(obj.subset) <- assay
    if (assay == "RNA") {
      obj.subset <- NormalizeData(obj.subset, verbose = FALSE)
      obj.subset <- ScaleData(obj.subset, features = genes, verbose = FALSE)
    } else if (assay == "SCT") {
      obj.subset <- GetResidual(obj.subset, features = genes)
    }
    p1 <- DoHeatmap(obj.subset, features = genes, slot = "scale.data") +
      theme(axis.text.y = element_text(size = 0)) +
      xlab(paste0("accuracy.1: ", round(accuracy.1, 2))) +
      NoLegend()
    p2 <- ggplot(test.data, aes(predscore.1, predscore.logistic)) +
      geom_point() +
      xlab("Pred.score.1") +
      ylab("Logistic.pred.score") +
      theme_bw()
    p3 <- ggplot(test.data, aes(predscore.2, predscore.logistic)) +
      geom_point() +
      xlab("Pred.score.2") +
      ylab("Logistic.pred.score") +
      theme_bw()
    layout <- "
  AB
  AC"
    p <- p1 + p2 + p3
    p <- p + patchwork::plot_layout(design = layout)

    return(list(
      test.data = test.data, accuracy.1 = accuracy.1,
      accuracy.2 = accuracy.2,
      plot = p
    ))
  }


  return(list(
    test.data = test.data, accuracy.1 = accuracy.1,
    accuracy.2 = accuracy.2
  ))
}

#' Compare accuracy of prediction of two different annotation methods
#'
#' This method uses consensus prediction between the methods as a ground truth to find markers for
#' each of the two celltypes where the two methods have a consensus and then uses these markers
#' to predict the most likely celltype of cells where the predictions differ
#' @param object Seurat object
#' @param celltype.1 Celltype 1
#' @param celltype.2 Celltype 2
#' @param predcol.1 Metadata column name for first method's predictions
#' @param predcol.2 Metadata column name for second methods' predictions
#' @param predcol.1 Metadata column name for first method's prediction scores
#' @param predcol.2 Metadata column name for second methods' prediction scores
#' @param min.diff Minimum number of differences to filter celltype differences on
#' @importFrom magrittr %>%
#' @importFrom dplyr filter arrange desc top_n bind_rows
#' @export
GetBestPrediction <- function(object,
                              predcol.1,
                              predcol.2,
                              scorecol.1,
                              scorecol.2,
                              min.diff = 10,
                              ...) {
  pred.diffs <- SummarisePredDiffs(object, predcol.1 = predcol.1, predcol.2 = predcol.2, min.diff.cells = min.diff, diff.only = T)
  # pred.diffs
  possible_combinations <- combinat::combn(unique(as.character(object@meta.data[, predcol.1])), 2)

  pred_diffs_filtered <- pred.diffs %>% filter(Freq >= min.diff)
  pred_diffs_filtered <- pred_diffs_filtered[pred_diffs_filtered[, predcol.1] %in% unique(possible_combinations[1, ]), ]
  accuracies <- list()
  plots <- list()
  test.data <- list()
  for (index in 1:nrow(pred_diffs_filtered)) {
    celltype.1 <- pred_diffs_filtered[index, 1]
    celltype.2 <- pred_diffs_filtered[index, 2]
    data <- ComparePredictions(
      object = object,
      celltype.1 = celltype.1,
      celltype.2 = celltype.2,
      predcol.1 = predcol.1,
      predcol.2 = predcol.2,
      scorecol.1 = scorecol.1,
      scorecol.2 = scorecol.2,
      ...
    )

    test.data[[index]] <- data[["test.data"]]
    accuracies[[index]] <- data.frame(
      acc.1 = data[["accuracy.1"]],
      acc.2 = data[["accuracy.2"]],
      celltype.1 = celltype.1,
      celltype.2 = celltype.2,
      differences = pred_diffs_filtered[index, 3]
    )
  }

  accuracy_df <- bind_rows(accuracies)
  accuracy_df$index <- rownames(accuracy_df)

  test_df <- bind_rows(test.data)
  test_df$index <- rownames(test_df)

  return(list(accuracy = accuracy_df, predictions = test_df))
}



#' Plot confusion matrices for prediction
#' @export
#' @importFrom ggplot2 geom_tile scale_fill_gradient theme_bw theme
PredictionComparisonHeatmap <- function(predictions.1, predictions.2) {
  common_cells <- intersect(rownames(predictions.1), rownames(predictions.2))

  df <- data.frame(predictions.1 = predictions.1[common_cells, 1], predictions.2 = predictions.2[common_cells, 1])

  predictions.compare <- table(df$predictions.1, df$predictions.2)
  predictions.compare <- predictions.compare / rowSums(predictions.compare) # normalize for number of cells in each cell type
  predictions.compare <- as.data.frame(predictions.compare)
  p <- ggplot(predictions.compare, aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    scale_fill_gradient(
      name = "Fraction of cells",
      low = "white", high = "blue"
    ) +
    ylab("predictions.1") +
    ylab("predictions.2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  return(p)
}
