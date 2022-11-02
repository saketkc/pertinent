#' importFrom GenomicRanges findOverlaps
SummarisePeakOverlap <- function(query, subject,
                                 query_name = "query", subject_name = "subject",
                                 query_peak_col = "peak", subject_peak_col = "peak") {
  hits <- findOverlaps(query = query, subject = subject)
  overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
  percentOverlap <- width(overlaps) / width(subject[subjectHits(hits)])

  query@elementMetadata[[paste0(subject_name, "_", subject_peak_col)]] <- NA
  query@elementMetadata[[paste0(subject_name, "_", subject_peak_col)]][queryHits(hits)] <- subject@elementMetadata[[subject_peak_col]][subjectHits(hits)]
  query@elementMetadata[[paste0("percent_overlap_with_", subject_name)]] <- NA
  query@elementMetadata[[paste0("percent_overlap_with_", subject_name)]][queryHits(hits)] <- percentOverlap

  subject@elementMetadata[[paste0(query_name, "_", subject_peak_col)]] <- NA
  subject@elementMetadata[[paste0(query_name, "_", subject_peak_col)]][subjectHits(hits)] <- query@elementMetadata[[query_peak_col]][queryHits(hits)]
  subject@elementMetadata[[paste0("percent_overlap_with_", query_name)]] <- NA
  subject@elementMetadata[[paste0("percent_overlap_with_", query_name)]][subjectHits(hits)] <- percentOverlap

  return(list(query = query, subject = subject))
}
