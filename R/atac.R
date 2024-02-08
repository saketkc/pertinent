#' import Signac
#' @export
AggFragObj <- function(fragpath, meta.data, genome) {
  ncounts <- CountFragments(fragments = fragpath)
  cells <- ncounts[ncounts$frequency_count > 500, "CB"]
  frags <- CreateFragmentObject(path = fragpath, cells = cells)

  agg_bins <- AggregateTiles(
    object = frags, genome = seqlengths(genome), min_counts = 5,
    cells = cells, binsize = 5000
  )
  col_sums <- colSums(agg_bins)
  row_sums <- rowSums(agg_bins)
  # agg_bins <- agg_bins[row_sums < 10000, col_sums < 5000]
  chrom_assay <- CreateAssayObject(counts = agg_bins, min.cells = 5, min.features = 200)
  obj <- CreateSeuratObject(counts = chrom_assay, assay = "tiles", )


  obj <- RunTFIDF(obj, scale.factor = median(obj$nCount_tiles))
  VariableFeatures(obj) <- rownames(obj)
  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = "lsi", dims = 2:25, verbose = F)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:25)
  obj <- FindClusters(obj, algorithm = 3)

  obj <- AddMetaData(object = obj, metadata = meta.data)

  return(obj)
}

#' @export
AggFragObj.regions <- function(fragpath, atac.features.granges, meta.data, frequency_count = 500,
                               binsize = 5000, dims = 2:25, min.features = 200, min.cells = 5, row_sums_filter = NULL) {
  atac <- FeatureMatrix(fragments = CreateFragmentObject(fragpath), features = atac.features.granges)

  chrom_assay <- CreateAssayObject(counts = atac, min.cells = min.cells, min.features = min.features)
  obj <- CreateSeuratObject(counts = chrom_assay, assay = "tiles", )


  obj <- RunTFIDF(obj, scale.factor = median(obj$nCount_tiles))
  # VariableFeatures(obj) <- rownames(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q0")

  obj <- RunSVD(obj)
  obj <- RunUMAP(obj, reduction = "lsi", dims = dims, verbose = F)
  obj <- FindNeighbors(obj, reduction = "lsi", dims = dims)
  obj <- FindClusters(obj, algorithm = 3)

  obj <- AddMetaData(object = obj, metadata = meta.data)

  return(obj)
}



#' Perform correction for Tn5 cut site
#' Increase forward mapping start positions by 4bp
#' Decrease reverse mapping end positions by 5bp.
#' @param df A bed red in as a dataframe
#' @export
DoTn5Correction <- function(df) {
  df.shifted <- df
  df.shifted$start[df$strand == "+"] <- df$start[df.shifted$strand == "+"] + 4
  df.shifted$end[df$strand == "-"] <- df$end[df.shifted$strand == "-"] - 5
  return(df.shifted)
}

#' Convert a given input bed pe file to a fragment file
#' to make it comptabile with downstream scATAC-seq workflows
#' @param bedpe.path Path to input bedpe file (can be gzipped)
#' @param fragment.path Path to output fragment file (gz compressed; sorted and tabix indexed)
#' @importFrom dplyr select filter arrange rowwise mutate
#' @importFrom data.table fread fwrite
#' @export
ConvertBedpeToFragment <- function(bedpe.path,
                                   fragment.path,
                                   chroms.to.keep = c(paste0("chr", seq(1, 22)), paste0("chr", c("X", "Y")))) {
  df <- fread(file = bedpe.path)
  colnames(df) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "barcode", "score", "strand1", "strand2")

  df.filtered <- df %>% filter(chrom1 == chrom2)

  df.long1 <- df.filtered %>% select(chrom1, start1, end1, barcode, score, strand1)
  colnames(df.long1) <- c("chrom", "start", "end", "barcode", "score", "strand")

  df.long2 <- df.filtered %>% select(chrom2, start2, end2, barcode, score, strand2)
  colnames(df.long2) <- c("chrom", "start", "end", "barcode", "score", "strand")

  df.filtered2 <- df %>%
    filter(chrom1 == chrom2) %>%
    filter(strand1 != strand2)

  # now separate them
  df.long1 <- df.filtered2 %>% select(chrom1, start1, end1, barcode, score, strand1)
  colnames(df.long1) <- c("chrom", "start", "end", "barcode", "score", "strand")

  df.long2 <- df.filtered2 %>% select(chrom2, start2, end2, barcode, score, strand2)
  colnames(df.long2) <- c("chrom", "start", "end", "barcode", "score", "strand")



  df.long1.corrected <- DoTn5Correction(df = df.long1)
  df.long2.corrected <- DoTn5Correction(df = df.long2)

  df.long1.corrected1 <- df.long1.corrected
  colnames(df.long1.corrected1) <- paste0(colnames(df.long1.corrected), "1")

  df.long2.corrected1 <- df.long2.corrected
  colnames(df.long2.corrected1) <- paste0(colnames(df.long2.corrected1), "2")

  df.long1long2 <- cbind(df.long1.corrected1, df.long2.corrected1) %>%
    rowwise() %>%
    mutate(
      chrom = chrom1,
      start = min(start1, start2, end1, end2),
      end = max(start1, start2, end1, end2),
      barcode = barcode1,
      score = score1
    )
  df.fragment <- df.long1long2 %>%
    select(chrom, start, end, barcode, score) %>%
    arrange(chrom, start, end, barcode, desc(score))

  if (!is.null(x = chroms.to.keep)) {
    df.fragment <- df.fragment %>% filter(chrom %in% chroms.to.keep)
  }
  df.fragment.gr <- GenomicRanges::sort(GenomicRanges::makeGRangesFromDataFrame(df = df.fragment, keep.extra.columns = T, starts.in.df.are.0based = T))
  df.fragment.sorted <- df.fragment.gr %>% as.data.frame() %>% select(seqnames, start, end, barcode, score)

  fwrite(df.fragment.sorted, file = fragment.path, col.names = FALSE, append = FALSE, sep = "\t")
  Rsamtools::bgzip(
    file = fragment.path,
    dest = sprintf("%s.bgz", sub("\\.gz$", "", fragment.path))
  )
  Rsamtools::indexTabix(file = sprintf("%s.bgz", sub("\\.gz$", "", fragment.path)), format = "bed")
}
