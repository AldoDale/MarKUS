#' Function to merge fastq files from different strands
#'
#' @param path Path of the directory with the fastq files
#' @param pattern Pattern to use to search for the files and to parse the sample names
#'
#' @return a data.frame with sample names and merged files path
#' @export
#'
#' @examples merge_fastq(path = ".", pattern = "*_nocont.fastq.1.gz")
setGeneric("merge_fastq", function(path, pattern)
  standardGeneric("merge_fastq"))

setMethod("merge_fastq",
          signature("character"), function(path = '.', pattern = '*_nocont.fastq.1.gz') {
            files <- Sys.glob(file.path(path, pattern))
            results <-
              data.frame(
                sample = character(),
                merged_path = character(),
                stringsAsFactors = FALSE
              )
            for (f1 in files) {
              base <- sub('\\.fastq\\.1\\.gz$', '', basename(f1))
              f2 <- file.path(path, paste0(base, '.fastq.2.gz'))
              out <- file.path(path, paste0(base, '_merged.fastq.gz'))
              cmd <- sprintf("bash -c %s",
                             shQuote(sprintf(
                               'zcat %s %s | gzip > %s',
                               shQuote(f1),
                               shQuote(f2),
                               shQuote(out)
                             )))
              message('Merging ',
                      basename(f1),
                      ' + ',
                      basename(f2),
                      ' to ',
                      basename(out))
              status <- system(cmd)
              if (status != 0)
                stop('merge_fastq failed on ', base)
              results <-
                rbind(results,
                      data.frame(
                        sample = base,
                        merged_path = out,
                        stringsAsFactors = FALSE
                      ))
            }
            return(results)
          })
