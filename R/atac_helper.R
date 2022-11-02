#' import Signac
#' @export
AggFragObj <- function(fragpath, meta.data, genome) {
  ncounts <- CountFragments(fragments = fragpath)
  cells <- ncounts[ncounts$frequency_count > 500, "CB"]
  frags <- CreateFragmentObject(path = fragpath, cells = cells)

  agg_bins <- AggregateTiles(
    object = frags, genome = seqlengths(genome), min_counts = 5,
    cells = cells, binsize = 5000
  )
  col_sums <- colSums(agg_bins)
  row_sums <- rowSums(agg_bins)
  # agg_bins <- agg_bins[row_sums < 10000, col_sums < 5000]
  chrom_assay <- CreateAssayObject(counts = agg_bins, min.cells = 5, min.features = 200)
  obj <- CreateSeuratObject(counts = chrom_assay, assay = "tiles", )


  obj <- RunTFIDF(obj, scale.factor = median(obj$nCount_tiles))
  VariableFeatures(obj) <- rownames(obj)
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = "lsi", dims = 2:25, verbose = F)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:25)
  obj <- FindClusters(obj, algorithm = 3)

  obj <- AddMetaData(object = obj, metadata = meta.data)

  return(obj)
}

#' @export
AggFragObj.regions <- function(fragpath, atac.features.granges, meta.data, frequency_count = 500,
                               binsize = 5000, dims = 2:25, min.features = 200, min.cells = 5, row_sums_filter = NULL) {
  atac <- FeatureMatrix(fragments = CreateFragmentObject(fragpath), features = atac.features.granges)

  chrom_assay <- CreateAssayObject(counts = atac, min.cells = min.cells, min.features = min.features)
  obj <- CreateSeuratObject(counts = chrom_assay, assay = "tiles", )


  obj <- RunTFIDF(obj, scale.factor = median(obj$nCount_tiles))
  # VariableFeatures(obj) <- rownames(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q0")

  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = "lsi", dims = dims, verbose = F)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, algorithm = 3)

  obj <- AddMetaData(object = obj, metadata = meta.data)

  return(obj)
}
