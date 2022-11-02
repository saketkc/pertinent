#' Convert a distance matrix to a data frame
#' @export
DistToDF <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B) - 1):1),
    value = as.vector(inDist)
  )
}

#' Create a pseudobulk of any seurat object with replicates
#' @importFrom Seurat AddMetaData AggregateExpression NormalizeData
#'
#' @export
CreateKPseudoBulk <- function(object, ident.name = "subclass.predicted.id", k_replicates = 3, seed = 42) {
  set.seed(seed)
  metadata <- object@meta.data
  pred_celltypes <- metadata[[ident.name]]
  celltype_composition <- table(pred_celltypes)
  replicate_labels <- paste0("replicate", seq(k_replicates))

  generated_labels <- lapply(names(celltype_composition),
    FUN = function(x) sample(replicate_labels, size = celltype_composition[[x]], replace = TRUE)
  )
  names(generated_labels) <- names(celltype_composition)

  metadata$replicate <- ""
  for (celltype in names(generated_labels)) {
    metadata_sub <- metadata[metadata[[ident.name]] == celltype, ]
    metadata[rownames(metadata_sub), "replicate"] <- generated_labels[[celltype]]
  }
  # for each celltype generate labels for each cell from replicate_labels such that
  # each celltype gets similar number of celltype labels

  composition <- as.data.frame(table(metadata[[ident.name]], metadata$replicate))
  colnames(composition) <- c("celltype", "replicate", "ncells")

  object <- AddMetaData(object, metadata = metadata)
  Idents(object) <- ident.name
  agg_obj <- AggregateExpression(object, assays = "RNA", slot = "counts", group.by = c(ident.name, "replicate"), return.seurat = T)
  metadata_celltypes <- stringr::str_split_fixed(colnames(agg_obj), "_", 2)[, 1]
  metadata_replicates <- stringr::str_split_fixed(colnames(agg_obj), "_", 2)[, 2]
  metadata <- data.frame(celltype = metadata_celltypes, replicate = metadata_replicates, sample_name = colnames(agg_obj))
  metadata <- as.data.frame(dplyr::left_join(metadata, composition, by = c("celltype", "replicate")))
  rownames(metadata) <- metadata$sample_name

  agg_obj <- AddMetaData(agg_obj, metadata = metadata)
  agg_obj <- NormalizeData(agg_obj)


  return(agg_obj)
}

#' Calculate average distance given a set of replicate distances
#' @export
AverageReplicates <- function(dist.matrix) {
  if (class(dist.matrix)[1] == "dist") {
    dists.df <- DistToDF(inDist = dist.matrix)
  } else {
    dists.df <- dist.matrix %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(-rowname)
    colnames(dists.df) <- c("row", "col", "value")
  }


  dists.df$row_replicate <- stringr::str_split_fixed(dists.df$row, pattern = "_", n = 3)[, 3]
  dists.df$col_replicate <- stringr::str_split_fixed(dists.df$col, pattern = "_", n = 3)[, 3]

  dists.df$row2 <- stringr::str_split_fixed(dists.df$row, pattern = "_replicate", n = 2)[, 1]
  dists.df$col2 <- stringr::str_split_fixed(dists.df$col, pattern = "_replicate", n = 2)[, 1]

  dists.df$row2 <- gsub("_", "-", dists.df$row2)
  dists.df$col2 <- gsub("_", "-", dists.df$col2)


  dists.df.mean <- dists.df %>%
    group_by(row2, col2) %>%
    summarise(value = mean(value))
  # return(dists.df.mean)

  dists.mean <- usedist::pivot_to_numeric_matrix(data = dists.df.mean, row2, col2, value)
  # dists.mean <- as.dist(dists.mean)
  return(dists.mean)
}
