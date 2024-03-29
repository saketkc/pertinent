#' Convert RDS to H5AD
#'
#' @param filepath Path to RDS file
#' @importFrom Seurat DietSeurat GetAssayData SetAssayData
#' @export
RDStoH5AD <- function(filepath) {
  h5seurat <- gsub(".rds", ".h5Seurat", filepath)
  query <- readRDS(filepath)
  DefaultAssay(query) <- "RNA"
  query <- DietSeurat(query, assays = "RNA", scale.data = FALSE, dimreducs = NULL, graphs = NULL)

  query <- SetAssayData(query, slot = "data", new.data = GetAssayData(query, slot = "counts", assay = "RNA"))
  SeuratDisk::SaveH5Seurat(query, filename = h5seurat, overwrite = TRUE, verbose = FALSE)
  SeuratDisk::Convert(source = h5seurat, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
}
