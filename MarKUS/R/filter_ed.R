#' Function to filter based on edit distances
#'
#' @param x a list of vectors
#' @param threshold minimum number of edits to reach the most similar sequence
#'
#' @return a named list
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom stringdist stringdistmatrix
#' @importFrom BiocParallel bplapply MulticoreParam
#' @export
#' @examples \dontrun{(x, threshold = 3, method = "lv", chunk_size  = 500, BPPARAM = MulticoreParam())}
filter_ed <- function(x, threshold = 3, method = "lv", chunk_size = 500,
                      BPPARAM = BiocParallel::MulticoreParam()) {

  if (!is.list(x) || any(!sapply(x, is.numeric))) {
    stop("x must be a named list of named numeric vectors.")
  }

  filt_final <- list()

  for (group in names(x)) {
    message("Processing ", group)

    kmer_vec <- names(x[[group]])
    this_set <- Biostrings::DNAStringSet(kmer_vec)
    N <- length(this_set)

    if (N == 0) {
      filt_final[[group]] <- character(0)
      next
    }

    idx_chunks <- if (chunk_size >= N) {
      list(seq_len(N))
    } else {
      split(seq_len(N), ceiling(seq_len(N) / chunk_size))
    }

    keep <- BiocParallel::bplapply(idx_chunks, function(chunk_idx) {
      chunk <- as.character(this_set[chunk_idx])
      rest <- as.character(this_set[-chunk_idx])
      mat <- stringdist::stringdistmatrix(chunk, rest, method = method)
      keep_idx <- apply(mat, 1, function(d) all(d > threshold))
      chunk[keep_idx]
    }, BPPARAM = BPPARAM)

    filt_final[[group]] <- unlist(keep, use.names = FALSE)
  }

  return(filt_final)
}
