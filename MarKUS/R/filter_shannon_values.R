#' Function to filter sequences based on their Shannon value
#'
#' @param x a named vector
#'
#' @return a named vector
#' @export
#'
#' @examples \dontrun{filter_shannon_values(x, threshold = 1.8)}
setGeneric("filter_shannon_values", function(x, threshold)
  standardGeneric("filter_shannon_values"))

setMethod("filter_shannon_values",
          signature("vector"), function(x, threshold) {
            filter_one <- function(v) {
              if (!is.numeric(v) || !is.atomic(v)) {
                stop("Each component must be a numeric vector")
              }
              if (is.null(names(v))) {
                stop("Each numeric vector must have names (the sequences)")
              }
              v[v > threshold]
            }

            # Dispatch on input type
            if (is.numeric(x) && is.atomic(x)) {
              # single vector
              return(list(filtered = filter_one(x)))
            }
            if (is.list(x) && !is.null(names(x))) {
              # list of vectors
              out <- lapply(x, filter_one)
              return(out)
            }

            stop("`x` must be either a named numeric vector or a named list of them")
          }
)

