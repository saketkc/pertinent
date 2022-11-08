
ReplaceMissingEntries <- function(df, col.source, col.target) {
  df[[col.target]][is.na(x = df[[col.target]])] <- df[[col.source]][is.na(x = df[[col.target]])]
  df[[col.target]][df[[col.target]] == ""] <- df[[col.source]][df[[col.target]] == ""]
  return(df)
}

#' Get longest isoform per gene
#' @importFrom GenomicFeatures  makeTxDbFromGFF transcriptsBy transcriptLengths
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols splitAsList
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange filter rename select
#' @export
ParseGTF <- function(gtf) {
  # gtf <- "/brahms/choudharys/data/Parsebio_Reference/Mus_musculus.GRCm39.108.gtf.gz"
  gtf.gr <- import(con = gtf, format = "gtf")

  gene_map <- mcols(x = gtf.gr)[, c("gene_id", "gene_name", "gene_biotype")] %>%
    as.data.frame() %>%
    unique() %>%
    arrange(gene_id)
  gene_map <- ReplaceMissingEntries(gene_map, col.source = "gene_id", col.target = "gene_name")

  tx_map <- mcols(x = gtf.gr)[, c("transcript_id", "gene_id", "gene_name", "gene_biotype")] %>%
    as.data.frame() %>%
    filter(!is.na(transcript_id)) %>%
    unique() %>%
    arrange(gene_id, transcript_id)
  tx_map <- ReplaceMissingEntries(tx_map, col.source = "gene_id", col.target = "gene_name")

  txdb <- suppressWarnings(expr = makeTxDbFromGFF(file = gtf, format = "gtf"))
  grl <- transcriptsBy(txdb, by = "gene")

  tx_lengths <- transcriptLengths(txdb = txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE) %>%
    rename(transcript_id = tx_name) %>%
    select(-tx_id)
  tx_lengths <- tx_lengths %>%
    left_join(y = tx_map) %>%
    select(transcript_id, gene_id, gene_name, gene_biotype, tx_len, cds_len, utr5_len, utr3_len) %>%
    arrange(gene_name, gene_id)

  return(list(tx_lengths = tx_lengths, txdb = txdb))
}
