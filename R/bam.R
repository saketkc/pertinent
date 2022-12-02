#' Split parsebio bam based on HexR or PolyT priming
#' @importFrom stringr str_split_fixed
#' @importFrom rtracklayer export
#' @importFrom Rsamtools ScanBamParam
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
#'
SplitParsebioBam <- function(file,
                             hexR.out = gsub(pattern = "\\.bam$",
                                             replacement = "_hexR.bam",
                                             x = file,
                                             ignore.case = TRUE),
                             polyT.out = gsub(pattern = "\\.bam$",
                                              replacement = "_polyT.bam",
                                              x = file,
                                              ignore.case = TRUE),
                             name.sep = "__",
                             readtype.num = 2,
                             barcode.tag = "CB",
                             cells.subset = NULL,
                             metadata = NULL,
                             group.by = NULL,
                             verbose = TRUE) {
  what <- c("flag", "mapq")

  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isSupplementaryAlignment = FALSE
  )
  param <- ScanBamParam(flag = flag, tag = c(barcode.tag), what = what)

  if (verbose) {
    message("Reading bam ...")
  }
  alignments <- readGAlignments(file = file, use.names = T, param = param)
  alignments <- alignments[alignments@elementMetadata$mapq >= 255]
  query_names <- alignments@NAMES
  if (verbose) {
    message("Extracting read names ...")
  }
  name_split <- str_split_fixed(string = query_names, pattern = name.sep, n = Inf)

  barcode <- mcols(x = alignments)[, barcode.tag]
  read_type <- name_split[, readtype.num]

  alignments@elementMetadata$barcode <- barcode
  alignments@elementMetadata$read_type <- read_type

  if (!is.null(x = cells.subset)) {
    if (verbose) {
      message("Filtering cells")
    }
    alignments <- alignments[alignments@elementMetadata$barcode %in% cells.subset]
  }

  if (verbose) {
    message("Filtering alignments")
  }
  hexR_alignments <- alignments[alignments@elementMetadata$read_type == "R"]
  polyT_alignments <- alignments[alignments@elementMetadata$read_type == "T"]

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
    message("Done.")
  }
  # df_summary <- df_filtered %>%
  #   group_by(GX, GN, read_type) %>%
  #   summarise(n_reads = n()) %>%
  #   arrange(GX, GN, read_type) %>%
  #   rename(gene_id = GX, gene_name = GN)
  return(df_filtered)
}

#' Split any bam
#' @importFrom dplyr filter left_join
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pbapply
#' @importFrom magrittr %>%
#' @importFrom rtracklayer export
#' @importFrom data.table fread setDF
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#' @importFrom S4Vectors mcols
#' @importFrom GenomicAlignments readGAlignments
#' @export
#'
SplitBam <- function(file, barcodes, out.dir, barcode.tag = "CB", verbose = TRUE) {
  flag <- scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE,
    isSupplementaryAlignment = FALSE
  )
  param <- ScanBamParam(tag = c(barcode.tag), flag = flag)
  if (verbose) {
    message("Reading bam ...")
  }
  alignments <- readGAlignments(file = file, use.names = T, param = param)

  bam_barcodes <- mcols(x = alignments)[, barcode.tag] %>% as.data.frame()
  colnames(x = bam_barcodes) <- "barcode"
  alignments@elementMetadata$barcode <- bam_barcodes$barcode

  given_barcodes <- setDF(x = fread(file = barcodes, header = FALSE))
  colnames(x = given_barcodes) <- c("barcode", "group")
  barcodes_to_write <- left_join(x = given_barcodes %>% distinct(), y = bam_barcodes %>% distinct(), by = "barcode")
  barcodes_to_write_list <- barcodes_to_write %>%
    group_by(group) %>%
    summarise(barcodes = paste0(barcode, collapse = ",")) %>%
    pull(barcodes, group) %>%
    as.list()
  split_barcodes <- function(x, split = ",") {
    return(unlist(x = strsplit(x = x, split = split)))
  }
  barcodes_to_write_list <- lapply(barcodes_to_write_list, FUN = split_barcodes)
  if (nbrOfWorkers() > 1) {
    mylapply <- future_lapply
  } else {
    mylapply <- pblapply
  }
  results <- mylapply(X = names(x = barcodes_to_write_list), FUN = function(group) {
    alignments_shortlist <- alignments[alignments@elementMetadata$barcode %in% barcodes_to_write_list[[group]]]
    suppressMessages(expr = suppressWarnings(expr = export(object = alignments_shortlist, con = file.path(out.dir, paste0(group, ".bam")), format = "bam")))
  })
}


#' Convert BAM to bigwig
#' @importFrom  rtracklayer export.bw
#' @importFrom GenomicAlignments coverage readGAlignmentPairs readGAlignments
#' @export
bam2bw <- function(file,
                   bw.path = gsub(pattern = "\\.bam$",
                                  replacement = ".bw",
                                  x = file, ignore.case = TRUE),
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
