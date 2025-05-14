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
  if (is.null(names(x))) {
    stop("`shannon_values` must be a named list")
  }
  lapply(x, function(vals) {
    vals[vals > threshold]
  })
          }
)

