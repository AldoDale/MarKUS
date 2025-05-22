#' Function to filter sequences based on the presence of repeated patterns
#'
#' @param x a list of vectors
#' @param mode used to select if the search should find a pattern (pattern) or any repeated sequence (any)
#' @param pattern the pattern to search for if mod = "pattern"
#' @param min_repeats the minimum number of repeats to search for
#' @param kmer_length used to choose the repeated sequence length if mode = "any"
#'
#' @return a list of vectors
#'
#' @export

setGeneric("filter_repeated_seqs", function(x,
                                             mode = c("pattern", "any"),
                                             pattern = NULL,
                                             min_repeats = 2,
                                             kmer_length = 3)
  standardGeneric("filter_repeated_seqs"))
setMethod("filter_repeated_seqs",
          signature("vector"),
          function(x,
                   mode        = c("pattern", "any"),
                   pattern     = NULL,
                   min_repeats = 2,
                   kmer_length = 3) {
            
            mode <- match.arg(mode)
            
            ## ---- define the two checkers ----
            has_repeat_vec <- function(seqs) {
              patt <- paste0("(?:", pattern, "){", min_repeats, ",}")
              grepl(patt, seqs, perl = TRUE)
            }
            
            has_any_kmer_repeat_vec <- function(seqs) {
              vapply(seqs, function(s) {
                L <- nchar(s)
                if (kmer_length > L) return(FALSE)
                kmers <- unique(substring(
                  s,
                  first = 1:(L - kmer_length + 1),
                  last  = kmer_length:(L)
                ))
                any(vapply(kmers, function(km) {
                  hits <- gregexpr(km, s, fixed = TRUE)[[1]]
                  (hits[1] != -1L) && (length(hits) >= min_repeats)
                }, logical(1)))
              }, logical(1))
            }
            
            checker <- if (mode == "pattern") {
              if (is.null(pattern) || !nzchar(pattern)) {
                stop("`pattern` must be a non-NULL, non-empty string when mode = 'pattern'")
              }
              has_repeat_vec
            } else {
              has_any_kmer_repeat_vec
            }
            
            ## ---- apply per group, *inverting* to keep the *non*-repeats ----
            res <- lapply(names(x), function(grp) {
              seqs <- x[[grp]]
              keep <- !checker(seqs)
              seqs[keep]
            })
            names(res) <- names(x)
            res
          })
