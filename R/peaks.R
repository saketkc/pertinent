#' @export
ReadBed <- function(path) {
  df <- setDF(x = fread(file = path))
  colnames(df) <- c("chrom", "start", "end", "name", "score", "strand")
  return(df)
}

#' @export
ReadNarrowPeak <- function(path, cols_to_keep = c("chrom", "start", "end", "score")) {
  peak_read <- read.table(path,
    col.names = c(
      "chrom",
      "start", "end", "name", "score", "strand", "fold_change",
      "neg_log10pvalue_summit", "neg_log10qvalue_summit",
      "relative_summit_position"
    )
  )
  if (!is.null(x = cols_to_keep)) peak_read <- peak_read[, cols_to_keep]
  return(peak_read)
}

#' @export
ReadSummit <- function(path, cols_to_keep = c("chrom", "start", "end", "score")) {
  df <- setDF(x = fread(file = path))
  colnames(df) <- c("chrom", "start", "end", "name", "score")
  if (!is.null(x = cols_to_keep)) df <- df[, cols_to_keep]
  gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  return(gr)
}

#' @export
GetChromSizes <- function(bsgenome, keep.stdchroms = TRUE, filter.mt = TRUE) {
  chrom.gr <- as(object = seqinfo(x = bsgenome), Class = "GRanges")
  if (keep.stdchroms) {
    chrom.gr <- keepStandardChromosomes(x = chrom.gr)
  }
  if (filter.mt) {
    chrom.gr <- chrom.gr[!seqnames(x = chrom.gr) %in% c("MT", "chrM", "chrMT")]
  }
  return(chrom.gr)
}
#' @export
FilterPeaks <- function(gr, filter.by = "score", decreasing = TRUE) {
  gr <- sort(x = sortSeqlevels(x = gr))
  # create a union of peaks
  gr.reduce <- reduce(gr, min.gapwidth = 0L, ignore.strand = TRUE)
  # find overlaps between given peak set and union of peaks
  overlaps <- findOverlaps(query = gr, subject = gr.reduce)
  # assing a cluster to given peaks based on the union of peaks
  mcols(gr)$cluster.reduce <- subjectHits(x = overlaps)

  if (!is.null(x = filter.by)) {
    # order peaks by score
    gr <- gr[order(mcols(x = gr)[, filter.by], decreasing = decreasing), ]
    # remove peaks that are assigned the same cluster retaining the best hit
    gr.dedup <- gr[!duplicated(x = mcols(x = gr)$cluster.reduce), ]
    # sort by chromosome and starts
    gr <- sort(x = sortSeqlevels(x = gr.dedup))
  } else {
    gr <- gr[!duplicated(mcols(x = gr)$cluster.reduce), ]
  }
  mcols(gr)$cluster.reduce <- NULL
  return(gr)
}

#' @export
IterativeFilterPeaks <- function(gr, filter.by = "score", decreasing = TRUE) {
  stopifnot(exprs = filter.by %in% colnames(x = mcols(x = gr)))
  i <- 0
  gr.initial <- gr # initial gr i
  while (length(gr.initial) > 0) {
    i <- i + 1
    # Keep strongest peaks
    gr.reduced <- FilterPeaks(gr = gr.initial, filter.by = filter.by, decreasing = decreasing)
    # remove called peaks
    gr.initial <- subsetByOverlaps(gr.initial, gr.reduced, invert = TRUE) # blacklist called cluster
    if (i == 1) {
      gr.all <- gr.reduced
    } else {
      gr.all <- c(gr.all, gr.reduced)
    }
  }
  return(gr.all)
}

#' Iteratively merge peaks extended around summits
#' Originally proposed in Corces and Granja 2018 (10.1126/science.aav1898)
#'
#' @export
#'
MergePeaksIterative <- function(summits.path.list = NULL,
                                summits.gr.list = NULL,
                                bsgenome = NULL,
                                blacklist = NULL, extend = 250,
                                min.spm = 1, ...) {
  if (is.null(x = summits.gr.list) & is.null(x = summits.path.list)) stop("Both summits.path.list and summits.gr.list cannot be null")
  chrom.sizes <- GetChromSizes(bsgenome = bsgenome, ...)

  if (!is.null(x = blacklist)) {
    blacklist <- import.bed(con = blacklist)
  } else {
    blacklist <- GRanges()
  }

  # 1. Read in summits
  if (is.null(summits.gr.list)) {
    summits.df <- lapply(X = summits.path.list, FUN = ReadSummit)
    summits.gr.list <- GRangesList(summits.df)
    remove(summits.df)
  }


  summits.gr.list <- lapply(X = seq_along(along.with = summits.gr.list), FUN = function(x) {
    summits.extend <- summits.gr.list[[x]] %>%
      # 2. Extend summits
      Extend(upstream = extend, downstream = extend, from.midpoint = TRUE) %>%
      # 3. Ensure they lie within chromosome boundaries
      subsetByOverlaps(., ranges = chrom.sizes, type = "within") %>%
      # 4. Remove blacklist overlaps
      subsetByOverlaps(., ranges = blacklist, type = "within", invert = TRUE) %>%
      # 5. Remove overlaps
      # IterativeFilterPeaks(., filter.by = "score", decreasing = TRUE)
      convergeClusterGRanges(., by = "score", decreasing = T)
    # 6. Add score per millionm
    mcols(x = summits.extend)$score <- mcols(x = summits.extend)$score / sum(mcols(x = summits.extend)$score) * 1e6
    return(summits.extend)
  })
  summits.gr.list <- GRangesList(summits.gr.list)
  # 7. Score by SPM and remove overlaps
  summits.gr.nonoverlapping <- summits.gr.list %>%
    unlist() %>%
    IterativeFilterPeaks(., filter.by = "score", decreasing = TRUE) %>%
    sortSeqlevels() %>%
    sort()

  # 8. Select based on SPM threshold
  summits.gr.nonoverlapping <- summits.gr.nonoverlapping[which(x = mcols(x = summits.gr.nonoverlapping)$score > min.spm), ]
  return(summits.gr.nonoverlapping)
}
