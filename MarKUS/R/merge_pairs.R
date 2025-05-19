#' Merge Paired-End FASTQ or FASTA files
#'
#' Merges paired-end sequencing files (R1 and R2) in various formats (.fastq, .fq, .fasta, .fa).
#'
#' @param path A character string specifying the directory containing the input files. Default is ".".
#' @param pattern_s1 Pattern to match first-strand (e.g., `*_1.fastq`, `*R1_001.fastq`) files. Default is "*_R1.fastq".
#' @param pattern_s2 Pattern to match second-strand (e.g., `*_2.fastq`, `*R2_001.fastq`) files. Default is "*_R2.fastq".
#' @param output_ext Desired file extension for the output.
#'
#' @return A data.frame with columns `sample` and `merged_path`.
#' @export
#'
#'@examples \dontrun{merge_pairs(path = "./path/to/files", pattern_s1 = "_R1.fastq", pattern_s2 = "_R2.fastq", output_ext = "fastq")}

setGeneric("merge_pairs", function(path, pattern_s1, pattern_s2, output_ext)
  standardGeneric("merge_pairs"))

setMethod("merge_pairs",
          signature("character"),
          function(path = ".",
                   pattern_s1 = "*_1.fastq.gz",
                   pattern_s2 = "*_2.fastq.gz",
                   output_ext = "fastq") {

            # Find R1 and R2 files
            files1 <- Sys.glob(file.path(path, pattern_s1))
            files2 <- Sys.glob(file.path(path, pattern_s2))

            if (length(files1) == 0L) stop("No files match pattern_s1 = '", pattern_s1, "'")
            if (length(files2) == 0L) stop("No files match pattern_s2 = '", pattern_s2, "'")

            suffix1 <- sub("^\\*", "", pattern_s1)
            suffix2 <- sub("^\\*", "", pattern_s2)

            results <- data.frame(
              sample = character(),
              merged_path = character(),
              stringsAsFactors = FALSE
            )

            for (f1 in files1) {
              base <- sub(paste0(suffix1, "$"), "", basename(f1))
              f2 <- file.path(path, paste0(base, suffix2))
              if (!file.exists(f2)) {
                stop("Mate file not found for '", base, "': expected '", f2, "'")
              }

              ext <- gsub("^\\.+", "", output_ext)
              out <- file.path(path, paste0(base, "_merged.", ext))

              compress <- grepl("\\.gz$", ext)

              cmd <- if (compress) {
                if (grepl("\\.gz$", f1)) {
                  sprintf("bash -c %s", shQuote(
                    sprintf("zcat %s %s | gzip > %s", shQuote(f1), shQuote(f2), shQuote(out))
                  ))
                } else {
                  sprintf("bash -c %s", shQuote(
                    sprintf("cat %s %s | gzip > %s", shQuote(f1), shQuote(f2), shQuote(out))
                  ))
                }
              } else {
                if (grepl("\\.gz$", f1)) {
                  sprintf("bash -c %s", shQuote(
                    sprintf("zcat %s %s > %s", shQuote(f1), shQuote(f2), shQuote(out))
                  ))
                } else {
                  sprintf("bash -c %s", shQuote(
                    sprintf("cat %s %s > %s", shQuote(f1), shQuote(f2), shQuote(out))
                  ))
                }
              }

              message("Merging ", basename(f1), " + ", basename(f2), " to ", basename(out))
              if (system(cmd) != 0) stop("Merge failed for ", base)

              results <- rbind(results, data.frame(
                sample = base,
                merged_path = out,
                stringsAsFactors = FALSE
              ))
            }

            return(results)
          }
)
