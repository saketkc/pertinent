#' Find the hamming distance between two strings
#'
#' @param input1 String
#' @param input2 String
#' @export
HammingDistance <- function(string1, string2) {
  s1 <- strsplit(x = string1, split = "")[[1]]
  s2 <- strsplit(x = string2, split = "")[[1]]
  if (length(x = s1) != length(x = s2)) stop("Inputs are not of the same length.")
  distance <- sum(s1 != s2)
  return(distance)
}

#' Return all possible 1bp mutated sequence
#' @param seq String
#' @return A vector of mutated sequences
#' @export
SingleMutations <- function(seq) {
  bases <- strsplit(x = "ACTG", split = "")[[1]]
  seq_chars <- strsplit(x = seq, split = "")[[1]]
  print(seq_chars)
  seq_length <- length(x = seq_chars)
  print(seq_length)
  mutated_sequences <- rep_len(x = "", length.out = (length(bases) - 1) * length(seq_chars))
  index <- 1
  for (i in 1:length(seq_chars)) {
    for (b in bases) {
      if (b != seq_chars[i]) {
        if (i == 1) {
          mutated_sequences[index] <- paste0(b, paste0(seq_chars[(i + 1):seq_length], collapse = ""), collapse = "")
        } else if (i == seq_length) {
          mutated_sequences[index] <- paste0(paste0(seq_chars[1:(i - 1)], collapse = ""), b, collapse = "")
        } else {
          mutated_sequences[index] <- paste0(paste0(seq_chars[1:(i - 1)], collapse = ""), b, paste0(seq_chars[(i + 1):seq_length], collapse = ""), collapse = "")
        }
        index <- index + 1
      }
    }
  }
  return(mutated_sequences)
}

AppendToDefaultDict <- function(obj.list, key, value) {
  if (key %in% names(x = obj.list)) {
    obj.list[[key]] <- c(obj.list[[key]], value)
  } else {
    obj.list[[key]] <- c(value)
  }
  return(obj.list)
}

#' Collpase kmer data for kmers that are within a hamming distance of 1.
#' This method is based on collapse_polyN function from Parsebio's splitpipe pipeline.
#' @param kmers A vector of kmers
#' @param counts A vecotr of counts associated with kmers
#' @description
#' @export
CollapseKmerCounts <- function(kmers, counts) {
  kmer_length <- length(x = strsplit(x = kmers[1], split = "")[[1]])
  # check if all kmers are of same length or not
  total_kmers <- length(kmers)
  kmer_seqs_dict <- rep_len(x = 0, length.out = total_kmers)
  names(kmer_seqs_dict) <- kmers

  ham_dict <- list()
  if (total_kmers > 50) {
    for (seq in kmers) {
      mutated_sequences <- SingleMutations(seq = seq)
      for (mut_seq in mutated_sequences) {
        if (!mut_seq %in% names(x = kmer_seqs_dict)) {
          ham_dict <- AppendToDefaultDict(ham_dict, mut_seq, seq)
        }
      }
    }
  } else {
    for (i in 1:(total_kmers - 1)) {
      for (j in (i + 1):total_kmers) {
        s1 <- kmers[i]
        s2 <- kmers[j]
        hamming_dist <- HammingDistance(string1 = s1, string2 = s2)
        if (hamming_dist <= 1) {
          ham_dict <- AppendToDefaultDict(ham_dict, s1, s2)
          ham_dict <- AppendToDefaultDict(ham_dict, s2, s1)
        }
      }
    }
  }
  kmer_counts <- counts
  names(x = kmer_counts) <- kmers
  kmer_counts <- as.list(x = kmer_counts)
  for (kmer in kmers) {
    cur_count <- kmer_counts[[kmer]]
    for (hm in ham_dict[[kmer]]) {
      # If there are kmers within hamming distance of 1 with more counts,
      # merge existing kmer into the one above and delete the current one
      if (hm %in% names(x = kmer_counts)) {
        if (kmer_counts[[hm]] > cur_count) {
          kmer_counts[[hm]] <- cur_count + kmer_counts[[hm]]
          kmer_counts[[kmer]] <- NULL
          break
        }
      }
    }
  }
  return(kmer_counts)
}
