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
MetaGeneCoverage <- function(bigwig.file,
                             gtf.data,
                             normalization.method = "max") {
  tx.pc <- gtf.data$tx_lengths %>%
    filter(gene_biotype == "protein_coding") %>%
    group_by(gene_id, gene_name) %>%
    top_n(wt = tx_len, n = 1) %>%
    filter(row_number() == 1)
  exons.gr <- exonsBy(x = gtf.data[["txdb"]], by = "tx", use.names = TRUE)[tx.pc$transcript_id]
  # bins <- matrix(data = NA, ncol = 100, nrow = length(exons.gr))
  # rownames(bins) <- names(x = exons.gr)
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- pblapply
  }
  results <- mylapply(X = exons.gr, FUN = function(transcript) {
    transcript <- sort(x = transcript)
    tx.strand <- unique(x = as.character(x = strand(x = transcript)))
    data <- unlist(x = import.bw(con = bigwig.file, which = transcript, as = "NumericList"))
    if (length(x = data) > 100) {
      if (tx.strand == "-") {
        data <- rev(x = data)
      }
      # Normalized in range of 0-1
      if (normalization.method == "max") {
        data_normalized <- data / max(data)
      } else if (normalization.method == "sum") {
        data_normalized <- data / sum(data)
      } else if (normalization.method == "median") {
        data_normalized <- data / median(data)
      } else if (normalization.method == "mean") {
        data_normalized <- data / mean(data)
      }

      df <- data.frame(x = 1:length(x = data_normalized), values = data_normalized)
      df$cut <- as.numeric(x = cut2(x = df$x, g = 100))
      df <- df %>%
        group_by(cut) %>%
        summarise(mean_va = mean(values), sum_va = sum(values)) %>%
        arrange(cut)
      if (nrow(df) >= 100) {
        return(df$mean_va)
      }
    }
  })
  bins <- do.call(what = rbind, args = results)
  bins_normalized <- colMeans(bins, na.rm = T)
  df_coverage <- data.frame(
    quantile = 1:length(bins_normalized),
    normalized_coverage = bins_normalized
  )
  return(df_coverage)
}

#' Pseudu bulk heatmap
#' @importFrom  pheatmap pheatmap
#' @importFrom Seurat AggregateExpression NormalizeData ScaleData
PseudoBulkHeatmap <- function(object, aggregate.by,
                              features.to.plot,
                              assays = "RNA",
                              cellheight = 20,
                              cellwidth = 20,
                              use.scaledata = FALSE,
                              idents.to.plot = NULL) {
  agg.counts <- AggregateExpression(object = object, assays = assays,
                                    group.by = aggregate.by, return.seurat = T)
  agg.counts <- NormalizeData(agg.counts, scale.factor = median(agg.counts$nCount_RNA))
  if (is.null(idents.to.plot)) {
    idents.to.plot <- colnames(agg.counts)
  }
  data <- agg.counts@assays$RNA@data[features.to.plot, idents.to.plot]
  scaled_data <- ScaleData(data)
  if (use.scaledata) {
    data.to.plot <- scaled_data
  } else {
    data.to.plot <- data
  }
  pheatmap(data.to.plot, cellheight = cellheight, cellwidth = cellwidth)
}
