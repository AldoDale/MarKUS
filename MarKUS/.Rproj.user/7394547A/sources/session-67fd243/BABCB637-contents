#' Function to extract and count k-mers
#'
#' @param obj data.frame with a "sample" column containing the sample names without file extention, and a "merged_path" column, with the full names of the merged files
#' @param kmer_size the chosen k-mer size
#' @param threads number of threads to use
#'
#' @return a data.frame with a "sample" column, a "count_db" column pointing at .jf files, and a "kmer_txt" column, pointing at the .txt files
#' @export
#'
#' @examples \dontrun{count_kmers(obj, kmer_size = 20, threads = 10)}
setGeneric("count_kmers", function(obj,
                                   kmer_size,
                                   threads)
  standardGeneric("count_kmers"))

setMethod("count_kmers",
          signature("data.frame"),
          function(obj,
                   kmer_size,
                   threads) {
            results <-
              data.frame(
                sample = character(),
                counts_db = character(),
                kmer_txt = character(),
                stringsAsFactors = FALSE
              )
            for (i in seq_len(nrow(obj))) {
              sample <- obj$sample[i]
              merged <- obj$merged_path[i]
              dir0 <- dirname(merged)
              db <- file.path(dir0, paste0(sample, '_counts_R.jf'))
              txt <- file.path(dir0, paste0(sample, '_kmer_counts_R.txt'))
              message('Counting k-mers for ', sample)

              inner_cmd <- sprintf(
                'jellyfish count -m %d -s 100M -t %d -C -o "%s" <(zcat "%s")',
                kmer_size,
                threads,
                db,
                merged
              )
              cmd1 <- sprintf('bash -o pipefail -c %s', shQuote(inner_cmd))
              status1 <- system(cmd1)
              if (status1 != 0)
                stop('jellyfish count failed on ', sample)

              cmd2 <- sprintf('jellyfish dump "%s" > "%s"',
                              db, txt)
              status2 <- system(cmd2)
              if (status2 != 0)
                stop('jellyfish dump failed on ', sample)
              message('Finished ', sample, '; output → ', basename(txt))
              results <-
                rbind(
                  results,
                  data.frame(
                    sample = sample,
                    counts_db = db,
                    kmer_txt = txt,
                    stringsAsFactors = FALSE
                  )
                )
            }
            return(results)
          }
)
