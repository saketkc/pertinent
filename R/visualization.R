#' Visualize bigwig coverage for a set of bigwigs over a gene
#' @importFrom magrittr %>%
#' @importFrom wiggleplotr plotCoverage getGenotypePalette
#' @importFrom dplyr select filter pull arrange desc
#' @importFrom GenomicFeatures exonsBy
VisualizeWigCoverage <- function(track.data, gene.name, gtf.data = NULL, gtf = NULL, exons.gr = NULL, tx.ids = NULL, ...) {
  if (is.null(x = gtf.data)) {
    # do some processing to get exons.gr
  }
  tx <- gtf.data[["tx"]]
  tx.lengths <- gtf.data[["tx_lengths"]]
  gene_of_interest <- tx.lengths %>%
    filter(gene_name == gene.name) %>%
    pull(gene_id) %>%
    unique()
  tx.gr <- tx[gene_of_interest]
  tx.gr.unlist <- unlist(x = tx.gr)
  tx_of_interest <- unique(x = tx.gr.unlist$tx_name)
  exons.gr <- exonsBy(x = gtf.data[["txdb"]], by = "tx", use.names = TRUE)[tx_of_interest]

  p <- suppressWarnings(expr = plotCoverage(exons.gr,
    track_data = track.data,
    heights = c(2, 1),
    fill_palette = getGenotypePalette(),
    ...
  ))
  return(p)
}
