
ReplaceMissingEntries <- function(df, col.source, col.target) {
  df[[col.target]][is.na(x = df[[col.target]])] <- df[[col.source]][is.na(x = df[[col.target]])]
  df[[col.target]][df[[col.target]] == ""] <- df[[col.source]][df[[col.target]] == ""]
  return(df)
}

#' Get longest isoform per gene
#' @importFrom GenomicFeatures  makeTxDbFromGFF transcriptsBy transcriptLengths
#' @importFrom rtracklayer import
#' @importFrom S4Vectors mcols splitAsList
#' @importFrom GenomeInfoDb
#' seqinfo seqinfo<-
#' seqnames seqnames<-
#' seqlevels seqlevels<-
#' sortSeqlevels
#' seqlengths seqlengths<-
#' isCircular isCircular<-
#' genome genome<-
#' @importFrom magrittr %>%
#' @importFrom BiocGenerics strand
#' @importFrom dplyr arrange distinct filter left_join rename select
#' @export
ParseGTF <- function(gtf, seqlevels_prefix = NULL) {
  gtf.gr <- import(con = gtf, format = "gtf")
  if (!is.null(x = seqlevels_prefix)) {
    seqlevels(x = gtf.gr) <- paste0(seqlevels_prefix, seqlevels(x = gtf.gr))
  }
  gene_map <- mcols(x = gtf.gr)[, c("gene_id", "gene_name", "gene_biotype")] %>% as.data.frame()
  gene_map$strand <- as.character(x = strand(x = gtf.gr))
  gene_map <- gene_map %>%
    distinct() %>%
    arrange(gene_id)
  gene_map$strand[gene_map$strand == "+"] <- 1
  gene_map$strand[gene_map$strand == "-"] <- -1
  gene_map <- ReplaceMissingEntries(gene_map, col.source = "gene_id", col.target = "gene_name")
  tx_map <- mcols(x = gtf.gr)[, c("transcript_id", "gene_id", "gene_name", "gene_biotype")] %>% as.data.frame()
  tx_map$strand <- as.character(x = strand(x = gtf.gr))
  tx_map <- tx_map %>%
    filter(!is.na(transcript_id)) %>%
    distinct() %>%
    arrange(gene_id, transcript_id)
  tx_map$strand[tx_map$strand == "+"] <- 1
  tx_map$strand[tx_map$strand == "-"] <- -1
  tx_map <- ReplaceMissingEntries(tx_map, col.source = "gene_id", col.target = "gene_name")
  txdb <- suppressWarnings(expr = makeTxDbFromGFF(file = gtf, format = "gtf"))
  if (!is.null(x = seqlevels_prefix)) {
    seqlevels(x = txdb) <- paste0(seqlevels_prefix, seqlevels(x = txdb))
  }
  grl <- transcriptsBy(txdb, by = "gene")
  tx_lengths <- transcriptLengths(txdb = txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE) %>%
    rename(transcript_id = tx_name) %>%
    select(-tx_id)
  tx_lengths <- tx_lengths %>%
    left_join(y = tx_map, by = c("transcript_id", "gene_id")) %>%
    select(transcript_id, gene_id, gene_name, gene_biotype, strand, tx_len, cds_len, utr5_len, utr3_len) %>%
    arrange(gene_name, gene_id)
  return(list(tx_lengths = tx_lengths, txdb = txdb, tx = grl))
}
