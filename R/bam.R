#' Split parsebio bam based on HexR or PolyT priming
#' @importFrom stringr str_split_fixed
#' @importFrom  rtracklayer export
#' @importFrom Rsamtools ScanBamParam
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
#'
SplitParsebioBam <- function(file,
                             hexR.out = gsub(pattern = "\\.bam$", replacement = "_hexR.bam", x = file, ignore.case = TRUE),
                             polyT.out = gsub(pattern = "\\.bam$", replacement = "_polyT.bam", x = file, ignore.case = TRUE),
                             name.sep = "__",
                             readtype.num = 2,
                             metadata = NULL,
                             group.by = NULL,
                             verbose = TRUE) {
  param <- ScanBamParam( # tag=c("nM",  "GX", "GN", "pN", "CB"),
    tag = c("CB")
  )
  if (verbose) {
    message("Reading bam ...")
  }
  alignments <- readGAlignments(file = file, use.names = T, param = param)

  query_names <- alignments@NAMES
  if (verbose) {
    message("Extracting read names ...")
  }
  name_split <- str_split_fixed(string = query_names, pattern = name.sep, n = Inf)

  barcode <- mcols(x = alignments)[, "CB"] # name_split[, barcode.num]
  read_type <- name_split[, readtype.num]

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
      if (verbose) {
        message("Creating bigwigs  ...")
      }
      export(object = hexR_alignments_subset, con = gsub(pattern = "\\.bam$", replacement = paste0("_", group, ".bam"), x = hexR.out, ignore.case = TRUE), format = "bam")
      export(object = polyT_alignmentss_subset, con = gsub(pattern = "\\.bam$", replacement = paste0("_", group, ".bam"), x = polyT.out, ignore.case = TRUE), format = "bam")
    }
  } else {
    if (verbose) {
      message("Creating bigwigs  ...")
    }
    export(object = hexR_alignments, con = hexR.out, format = "bam")
    export(object = polyT_alignments, con = polyT.out, format = "bam")
  }
}

#' Count reads in a parsebio output bam
#' @importFrom dplyr filter group_by summarise n rename arrange
#' @importFrom magrittr %>%
#' @importFrom stringr str_split
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
#'
CountReadsinParsebioBam <- function(file,
                                    name.sep = "__",
                                    readtype.num = 2,
                                    cells = NULL, verbose = TRUE) {
  what <- c("qname", "flag", "mapq")
  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isSupplementaryAlignment = FALSE
  )
  param <- ScanBamParam(
    flag = flag,
    tag = c("GX", "GN", "pN", "RE", "CB"),
    what = what
  )
  if (verbose) {
    message("Reading bam ...")
  }
  alignments <- readGAlignments(file = file, param = param)
  if (verbose) {
    message("Extracting read names ...")
  }
  df <- mcols(x = alignments) %>% as.data.frame()
  query_names <- df$qname
  name_split <- str_split_fixed(string = query_names, pattern = name.sep, n = Inf)
  read_type <- name_split[, readtype.num]
  if (verbose) {
    message("Creating mapping summary ...")
  }
  df$read_type <- read_type
  df$GN[df$GN == ""] <- df$GX[df$GN == ""]
  df_filtered <- df[(df$GX != "") | (df$GN != ""), ]

  if (!is.null(x = cells)) {
    df_filtered <- df_filtered %>% filter(CB %in% cells)
  }
  if (verbose) {
    message("Summarising ...")
  }
  # df_summary <- df_filtered %>%
  #   group_by(GX, GN, read_type) %>%
  #   summarise(n_reads = n()) %>%
  #   arrange(GX, GN, read_type) %>%
  #   rename(gene_id = GX, gene_name = GN)
  return(df_filtered)
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
    if (!is.null(x = normalization.method)) {
      if (tolower(x = normalization.method) == "cpm") {
        scale.factor <- as.numeric(x = 1000000 / n_reads)
        cov <- round(x = cov * scale.factor, digits = 3)
      } else if (tolower(x = normalization.method) == "rpkm") {
        # pass for now
        cov <- cov
      }
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
