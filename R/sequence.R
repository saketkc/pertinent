#' Fetch sequences from a DNAStringSet object
#' @param genomeseq A DNAStringSet object 
#' @param gr A granges object
#' @importFrom magrittr %>%
#' @importFrom Biostrings reverseComplement
#' @export
GetSequencesFromGenome <- function(genomeseq, gr){
  granges.df <- gr %>% as.data.frame()
  getsequences <- function(genomeseq, chr, start, end, strand = "+") {
    start <- as.numeric(start)
    end <- as.numeric(end)
    if (strand %in% c("+", "*")) {
      seq.char <- strsplit(x = as.character(genomeseq[[chr]]), split = "")[[1]]
    } else if (strand == "-") {
      seq.char <- strsplit(x = as.character(reverseComplement(genomeseq[[chr]])), split="")[[1]]
    }
    return(paste0(seq.char[start:end], collapse=""))
  }
  sequences <- apply(granges.df, MARGIN = 1, FUN = function(row) getsequences(genomeseq, row[["seqnames"]], row[["start"]], row[["end"]], row[["strand"]]))
  names(sequences) <- paste0(granges.df$seqnames, "-", granges.df$start, "-", granges.df$end, "-", granges.df$strand)
  return(sequences)
}

