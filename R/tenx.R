#' Convert a given integer to 2-bit and eventually ATCG..
#' sequence
#' Based on https://kb.10xgenomics.com/hc/en-us/articles/360004105372-How-do-you-decompress-the-2-bit-barcode-sequences-in-molecule-info-h5-file-
#'
#' @param seq.bit An integer (to be converted to 2-bit)
#' @return A character containing the leetters attro
#' @export
TwoBitToSeq <- function(seq.bit, length = 10) {
  # pad
  seq.bit.padded <- sprintf(paste0("%0", 2 * length, "d"), seq.bit)
  # convert to binary
  seq.bit.binary <- rev(x = as.numeric(x = intToBits(x = seq.bit.padded)))
  seq.bit.binary <- paste0(seq.bit.binary[-(1:(length(seq.bit.binary) - 2 * length))], collapse = "")
  nucleotide_list <- substring(text = seq.bit.binary, first = seq(1, nchar(seq.bit.binary), 2), last = seq(2, nchar(seq.bit.binary), 2))
  conversion_list <- list("00" = "A", "01" = "C", "10" = "G", "11" = "T")
  return(paste0(unlist(conversion_list[nucleotide_list]), collapse = ""))
}

#' Read 10X molecular_info.h5 file
#' @param filename Path to molecular_info.h5
#' @param version.10X "v2" or "v3"
#' @param make.umi.seq Whether to convert the integer umi values to a sequence (post 2-bit encoding/decoding)
#' @param return.raw Whether a raw value or just the read count summary at the cell barcode level
#' @export
Read10X_molinfo <- function(filename, version.10x = NULL, make.umi.seq = FALSE, return.raw = FALSE) {
  h5 <- hdf5r::H5File$new(filename = filename, mode = "r")
  h5.slots <- names(x = h5)
  if (is.null(x = version.10x)) {
    if ("barcode_idx" %in% h5.slots) {
      version.10x <- "3"
    } else {
      version.10x <- "2"
    }
  }

  if (version.10x == "3") {
    barcodes.all <- h5[["barcodes"]][]
    barcodes.idx <- h5[["barcode_idx"]][]
    # this is zero baed
    barcodes.name <- barcodes.all[barcodes.idx + 1L]
    readcounts <- h5[["count"]][]
    features.idx <- h5[["feature_idx"]][]
    features.lookup <- h5[["features/id"]][]
    # this is zero baed
    features.name <- features.lookup[features.idx + 1L]
  } else {
    # this is slow
    barcode.length <- 16 # for v2
    # barcode.length <- 14 # for v1
    barcodes.all <- h5[["barcode"]][]
    barcode.coorrected_reads <- h5[["barcode_corrected_reads"]][]

    barcodes.name <- rep("", length(barcodes.all))
    barcodes.name <- sapply(barcodes.all, TwoBitToSeq, length = barcode.length)
    # l <- lineprof(barcodes.name <- sapply(barcodes.all[1:100], TwoBitToSeq, length=barcode.length))
    # l
    readcounts <- h5[["reads"]][]
    features.idx <- h5[["gene"]][]

    features.lookup <- h5[["gene_ids"]][]
    gene_names <- h5[["gene_names"]][]
    features.name <- features.lookup[features.idx + 1L]
  }


  umi_int <- h5[["umi"]][]
  gem_group <- h5[["gem_group"]][]


  df <- data.frame(
    barcode = barcodes.name,
    gem_group = gem_group,
    feature = features.name,
    read_count = readcounts,
    umi_int = umi_int
  )
  df$barcode <- factor(x = df$barcode)
  total_reads <- tapply(df$read_count, df$barcode, sum)
  df.cells <- data.frame(total_reads)

  if (make.umi.seq == TRUE) {
    if (version.10x == "3") {
      # See: https://kb.10xgenomics.com/hc/en-us/articles/360004105372-How-do-you-decompress-the-2-bit-barcode-sequences-in-molecule-info-h5-file-
      umi.length <- 12
    } else {
      umi.length <- 10
    }
    df$umi_seq <- sapply(df$umi_int, TwoBitToSeq, length = umi.length)
  }

  if (return.raw) {
    return(list(molinfo = df, readinfo = df.cells))
  } else {
    return(df.cells)
  }
}
