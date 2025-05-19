#' @include classes.R
NULL

#' Function to calculate shannon entropy of sequences
#'
#' @param x named vector obtained from "get_shared_kmers" function
#'
#' @return a list containing a list of named vectors and a list of plots of shannon values distribution
#' @export
#'
#' @examples \dontrun{calculate_entropy(x)}
setGeneric("calculate_shannon", function(x)
  standardGeneric("calculate_shannon"))

setMethod("calculate_shannon",
          signature("dfOrVec"), function(x) {
            if (is.character(x) && !is.null(names(x))) {
              df <- data.frame(
                group = names(x),
                shared_path = unname(x),
                stringsAsFactors = FALSE
              )
            } else if (is.data.frame(x) &&
                       all(c("group","shared_path") %in% names(x))) {
              df <- x
            } else {
              stop("`x` must be either a named character vector or a data.frame with columns 'group' & 'shared_path'")
            }

            calc_ent <- function(seq) {
              bases <- strsplit(seq, split = "")[[1]]
              freqs <- table(bases) / length(bases)
              -sum(freqs * log2(freqs))
            }

            shannon_values <- list()
            shannon_plots <- list()

            for (i in seq_len(nrow(df))) {
              grp <- df$group[i]
              path <- df$shared_path[i]

              seqs <- readLines(path, warn = FALSE)
              if (length(seqs) == 0L) next

              ents <- vapply(seqs, calc_ent, numeric(1))
              shannon_values[[grp]] <- ents

              hist_df <- data.frame(shannon = ents, stringsAsFactors = FALSE)
              p <- ggplot2::ggplot(hist_df, ggplot2::aes(x = shannon)) +
                ggplot2::geom_histogram(binwidth = 0.05, fill = "grey70", color = "black") +
                ggplot2::labs(
                  title = paste("Shannon Entropy Distribution \u2013", grp),
                  x = "Shannon entropy",
                  y = "Count"
                ) +
                ggplot2::theme_minimal()
              shannon_plots[[grp]] <- p
            }

            list(
              shannon_values = shannon_values,
              shannon_plots = shannon_plots
            )
          }
          )
