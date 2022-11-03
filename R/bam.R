#' Split parsebio bam based on HexR or PolyT priming
#' @importFrom stringr str_split_fixed
#' @importFrom  rtracklayer export
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
#'
SplitParsebioBam <- function(file,
                             hexR.out = gsub(pattern = "\\.bam$", replacement = "_hexR.bam", x = file, ignore.case = TRUE),
                             polyT.out = gsub(pattern = "\\.bam$", replacement = "_polyT.bam", x = file, ignore.case = TRUE),
                             metadata = NULL,
                             group.by = NULL) {
  alignments <- readGAlignments(file = file, use.names = T)

  query_names <- alignments@NAMES
  name_split <- str_split_fixed(query_names, pattern = "__", n = 3)

  barcode <- name_split[, 1]
  read_type <- name_split[, 2]

  alignments@elementMetadata$barcode <- barcode
  alignments@elementMetadata$read_type <- read_type

  hexR_alignments <- alignments[read_type == "R"]
  polyT_alignments <- alignments[read_type == "T"]

  if (!is.null(metadata) && !is.null(group.by)) {
    groups <- unique(x = metadata[[group.by]])
    all_barcodes <- rownames(x = metadata)
    groupwise_barcode <- sapply(X = groups, FUN = function(x) all_barcodes[metadata[[group.by]] == x])
    for (group in groups) {
      hexR_alignments_subset <- hexR_alignments[hexR_alignments@elementMetadata$barcode %in% groupwise_barcode[[group]]]
      polyT_alignmentss_subset <- polyT_alignments[polyT_alignments@elementMetadata$barcode %in% groupwise_barcode[[group]]]
      export(object = hexR_alignments_subset, con = gsub(pattern = "\\.bam$", replacement = paste0("_", group, ".bam"), x = hexR.out, ignore.case = TRUE), format = "bam")
      export(object = polyT_alignmentss_subset, con = gsub(pattern = "\\.bam$", replacement = paste0("_", group, ".bam"), x = polyT.out, ignore.case = TRUE), format = "bam")
    }
  } else {
    export(object = hexR_alignments, con = hexR.out, format = "bam")
    export(object = polyT_alignments, con = polyT.out, format = "bam")
  }
}

#' Convert BAM to bigwig
#' @importFrom  rtracklayer export.bw
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
bam2bw <- function(file,
                   bw.path = gsub(pattern = "\\.bam$", replacement = ".bw", x = file, ignore.case = TRUE),
                   paired = FALSE,
                   strand = "both",
                   method = c("native", "deeptools"),
                   bamcoverage.path = NULL,
                   normalization.method = NULL,
                   cores = 1L) {
  method <- match.arg(method)
  if (method == "native") {
    alignments <- readGAlignments(file = file)
    cov <- coverage(x = alignments)
    n_reads <- sum(sum(cov))
    if (tolower(x = normalization.method) == "cpm") {
      scale.factor <- as.numeric(x = 1000000 / n_reads)
      cov <- round(x = cov * scale.factor, digits = 3)
    } else if (tolower(x = normalization.method) == "rpkm") {
      # not sure if this is correct
      cov <- cov
    }
    export.bw(object = cov, con = bw.path)
  } else if (method == "deeptools") {
    cmd <- paste(
      bamcoverage.path, "-b", file, "-o", bw.path,
      "--binSize 50 -of bigwig", "-p", cores
    )
    system(command = cmd)
  }
}
