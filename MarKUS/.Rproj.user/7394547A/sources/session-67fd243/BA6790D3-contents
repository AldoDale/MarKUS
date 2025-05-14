#' Function to filter based on edit distances
#'
#' @param x a list of vectors
#' @param threshold minimum number of edits to reach the most similar sequence
#'
#' @return a named list
#' @export
#'
#' @examples \dontrun{(x, threshold = 3, method = "lv", chunk_size  = 500, BPPARAM = MulticoreParam())}
setGeneric("filter_ed", function(x,
                                 threshold = 3,
                                 method = "lv",
                                 chunk_size = 500,
                                 BPPARAM = MulticoreParam())
  standardGeneric("filter_ed"))

setMethod("filter_ed",
          signature("vector"), function(x,
                                        threshold   = 3,
                                        method      = "lv",
                                        chunk_size  = 500,
                                        BPPARAM     = MulticoreParam()) {
            ## 1) Sanity‐checks
            if (is.null(names(x)) || any(names(x) == "")) {
              stop("`x` must be a named list of DNAStringSet objects")
            }
            if (!all(vapply(x, function(x)
              is(x, "DNAStringSet"), logical(1)))) {
              stop("All elements of `x` must be DNAStringSet objects")
            }

            result_list <- list()

            ## 2) Loop over each focal group
            for (grp in names(x)) {
              this_set   <- x[[grp]]
              other_sets <- DNAStringSet()
              for (o in setdiff(names(x), grp)) {
                other_sets <- c(other_sets, x[[o]])
              }
              other_chr <- as.character(other_sets)

    N <- length(this_set)
    ## split indices into chunks
    idx_chunks <- split(seq_len(N), ceiling(seq_len(N) / chunk_size))

    ## 3) Parallel chunk processing
    keep_list <- bplapply(idx_chunks, function(idxs) {
      # each worker sees `this_set`, `other_chr`, etc.
      chunk_chr <- as.character(this_set[idxs])
      Dmat      <- stringdistmatrix(chunk_chr,
                                    other_chr,
                                    method  = method,
                                    nthread = 1)
      # keep rows where all distances > threshold
      apply(Dmat, 1, function(r) all(r > threshold))
    }, BPPARAM = BPPARAM)

    ## combine logicals in original order
    keep_idx <- unlist(keep_list, use.names = FALSE)

    ## 4) collect survivors
    if (any(keep_idx)) {
      result_list[[grp]] <- DataFrame(
        sequence = this_set[keep_idx],
        group    = rep(grp, sum(keep_idx))
      )
    }
  }

  ## 5) bind and return
  if (length(result_list) == 0L) {
    return(DataFrame(sequence = DNAStringSet(), group = character()))
  }
  do.call(rbind, result_list)
          }
)
