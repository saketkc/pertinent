#' Get correlation per gene given two matrices (genes as rows, celltypes as columns)
#' @param mat1 Normalized matrix with genes as rows, cells as columns
#' @param mat2 Normalized matrix with genes as rows, cells as columns
#' @return A named list of genewise correlations.
#' @export
PerGeneCorrelation <- function(mat1, mat2) {
  if (dim(mat1)[2] != dim(mat2)[2]) {
    stop("Unequal column lengths")
  }
  if (dim(mat1)[1] != dim(mat2)[1]) {
    stop("Unequal gene lengths")
  }

  common_genes <- intersect(rownames(mat1), rownames(mat2))

  # sort to ensure columns are the same
  mat1 <- mat1[common_genes, sort(colnames(mat1))]
  mat2 <- mat2[common_genes, sort(colnames(mat2))]

  # we assume the columns are

  centered.mat1 <- mat1 - rowMeans(mat1)
  centered.mat2 <- mat2 - rowMeans(mat2)

  sd.mat1 <- sqrt(rowMeans(centered.mat1^2))
  sd.mat2 <- sqrt(rowMeans(centered.mat2^2))

  corr.rows <- rowMeans((centered.mat1 * centered.mat2) / (sd.mat1 * sd.mat2))
  return(corr.rows)
}
