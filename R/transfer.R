#' @importFrom Seurat FindTransferAnchors TransferData
#
#' @export
V3Transfer <- function(query, reference, refdata, reference.assay = "SCT", query.assay = "SCT",
                       normalization.method = "SCT", dims = 1:50) {
  transfer.anchors <- FindTransferAnchors(
    reference = reference, query = query,
    dims = dims, reduction = "cca", reference.assay = reference.assay, query.assay = query.assay,
    normalization.method = normalization.method, features = intersect(
      rownames(reference),
      rownames(query)
    ), n.trees = 20
  )

  predicted.labels <- TransferData(
    anchorset = transfer.anchors, refdata = refdata,
    weight.reduction = query[["pca"]], dims = dims
  )
  predicted.labels$predicted.id <- factor(predicted.labels$predicted.id, levels = sort(unique(as.character(predicted.labels$predicted.id))))

  return(list(transfer.anchors = transfer.anchors, predicted.labels = predicted.labels))
}

#' @import magrittr
#' @import Seurat
#
#' @export
SCTPrep <- function(rna_seu, dims = 1:50) {
  rna_seu <- SCTransform(
    object = rna_seu, assay = "RNA", method = "glmGamPoi",
    ncells = 5000, do.correct.umi = FALSE, do.scale = FALSE, do.center = TRUE
  ) %>%
    RunPCA(verbose = FALSE, npcs = max(dims)) %>%
    RunUMAP(dims = dims, verbose = FALSE) %>%
    FindNeighbors(dims = dims, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)
  return(rna_seu)
}



# Find anchors between query and reference
#' @import Seurat
#
#' @export
AziMuthMapping <- function(query, reference) {
  anchors <- FindTransferAnchors(
    reference = reference$map, query = query, k.filter = NA,
    reference.neighbors = "spca.annoy.neighbors", reference.assay = "SCT", query.assay = "SCT",
    reference.reduction = "spca", normalization.method = "SCT", features = intersect(
      rownames(x = reference$map),
      VariableFeatures(object = query)
    ), dims = 1:50, n.trees = 20, mapping.score.k = 100
  )

  query <- TransferData(
    reference = reference$map, query = query, dims = 1:50,
    anchorset = anchors, refdata = list(id = Idents(reference$map), impADT = GetAssayData(
      object = reference$map[["ADT"]],
      slot = "data"
    )), store.weights = TRUE
  )

  # Calculate the embeddings of the query data on the reference SPCA
  query <- IntegrateEmbeddings(
    anchorset = anchors, reference = reference$map,
    query = query, reductions = "pcaproject", reuse.weights.matrix = TRUE
  )

  # Calculate the query neighbors in the reference with respect to the integrated
  # embeddings
  query[["query_ref.nn"]] <- FindNeighbors(
    object = Embeddings(reference$map[["spca"]]),
    query = Embeddings(query[["integrated_dr"]]), return.neighbor = TRUE, l2.norm = TRUE
  )

  # The reference used in the app is downsampled compared to the reference on which
  # the UMAP model was computed. This step, using the helper function NNTransform,
  # corrects the Neighbors to account for the downsampling.
  query <- NNTransform(object = query, meta.data = reference$map[[]])

  # Project the query to the reference UMAP.
  query[["proj.umap"]] <- RunUMAP(
    object = query[["query_ref.nn"]], reduction.model = reference$map[["jumap"]],
    reduction.key = "UMAP_"
  )


  # Calculate mapping score and add to metadata
  query <- AddMetaData(
    object = query, metadata = MappingScore(anchors = anchors),
    col.name = "mapping.score"
  )
  return(query)
}
