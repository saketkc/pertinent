#' Convert RDS to H5AD
#'
#' @param filepath Path to RDS file
#' @importFrom Seurat DietSeurat GetAssayData SetAssayData
#' @importFrom  SeuratDisk SaveH5Seurat Convert
#'
#' @export
RDStoH5AD <- function(filepath) {
  h5seurat <- gsub(".rds", ".h5Seurat", filepath)
  DefaultAssay(query) <- "RNA"
  query <- DietSeurat(query, assays = "RNA", scale.data = FALSE, dimreducs = NULL, graphs = NULL)

  query <- SetAssayData(query, slot = "data", new.data = GetAssayData(query, slot = "counts", assay = "RNA"))
  SaveH5Seurat(query, filename = h5seurat, overwrite = TRUE, verbose = FALSE)
  Convert(source = h5seurat, dest = "h5ad", overwrite = TRUE, verbose = FALSE)
}
