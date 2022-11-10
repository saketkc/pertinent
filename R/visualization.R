#' Visualize bigwig coverage for a set of bigwigs over a gene
#' @importFrom BiocGenerics strand
#' @importFrom magrittr %>%
#' @importFrom wiggleplotr plotCoverage getGenotypePalette
#' @importFrom dplyr select filter pull arrange desc
#' @importFrom GenomicFeatures exonsBy
#' @export
VisualizeWigCoverage <- function(track.data, gene.name, gtf.data = NULL, gtf = NULL,
                                 exons.gr = NULL, tx.ids = NULL, heights = c(2, 1), ...) {
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

  p <- suppressWarnings(expr = plotCoverage(
    exons = exons.gr, transcript_annotations = tx.lengths,
    track_data = track.data,
    heights = heights,
    fill_palette = getGenotypePalette(),
    ...
  ))
  return(p)
}

#' Plot bigwig coverage as metagene
#' @importFrom BiocGenerics strand sort
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
#' @importFrom magrittr %>%
#' @importFrom rtracklayer import.bw
#' @importFrom dplyr arrange filter group_by top_n row_number summarise
#' @importFrom GenomicFeatures exonsBy
#' @importFrom Hmisc cut2
#' @export
MetaGeneCoverage <- function(bigwig.file, gtf.data) {
  tx.pc <- gtf.data$tx_lengths %>%
    filter(gene_biotype == "protein_coding") %>%
    group_by(gene_id, gene_name) %>%
    top_n(wt = tx_len, n = 1) %>%
    filter(row_number() == 1)
  exons.gr <- exonsBy(x = gtf.data[["txdb"]], by = "tx", use.names = TRUE)[tx.pc$transcript_id]
  bins <- matrix(data = NA, ncol = 100, nrow = length(exons.gr))
  rownames(bins) <- names(x = exons.gr)

  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- pblapply
  }
  results <- mylapply(X = names(x = exons.gr), FUN = function(txid) {
    transcript <- sort(x = exons.gr[[txid]])
    tx.strand <- unique(x = as.character(x = strand(x = transcript)))
    data <- unlist(x = import.bw(con = bigwig.file, which = transcript, as = "NumericList"))

    if (length(x = data)<100){
      next
    }
    if (tx.strand == "-") {
      data <- rev(x = data)
    }
    # Normalized in range of 0-1
    data_normalized <- data / max(data)
    df <- data.frame(x = 1:length(x=data_normalized), values = data_normalized)
    df$cut <- as.numeric(x = cut2(x = df$x, g = 100))
    df <- df %>%
      group_by(cut) %>%
      summarise(mean_va = mean(values), sum_va = sum(values)) %>%
      arrange(cut)
    if (nrow(df) >= 100) {
      bins[txid, ] <- df$mean_va
    }
  })

  bins_normalized <- colMeans(bins, na.rm = T)
  df_coverage <- data.frame(
    quantile = 1:length(bins_normalized),
    normalized_coverage = bins_normalized
  )
  return(df_coverage)
}
