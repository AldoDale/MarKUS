#' Function to extract unique k-mers
#'
#' @param obj a data.frame
#' @param data a metadata object with a "sample" column
#' @param group a column of the metadata to use as grouping factor
#' @param min_occurrence minimum number of samples in which k-mers must be found
#'
#' @return a data.frame
#' @export
#'
#' @examples \dontrun{get_shared_kmers(x = obj, data = meta, group= "treatment", min_occurrence = 4)}
setGeneric("get_shared_kmers", function(obj,
                                        data,
                                        group,
                                        min_occurrence)
  standardGeneric("get_shared_kmers"))

setMethod("get_shared_kmers",
          signature("data.frame"),
          function(obj, data, group, min_occurrence = 4) {
            mergeddf <- merge(obj, data, by = "sample")
            grsplit  <- split(mergeddf, mergeddf[[group]])

            out_files <- setNames(character(length(grsplit)), names(grsplit))

            for (grp in names(grsplit)) {
              dfg   <- grsplit[[grp]]
              files <- dfg$kmer_txt
              dir0  <- dirname(files[1])
              shared_out <-
                file.path(dir0, paste0(grp, "_shared_kmers_R.txt"))

              message("Working on group: ", grp)
              message("  Files: ", paste(basename(files), collapse = ", "))

              # 1) Extract unique kmers per file
              cmd1 <- sprintf("bash -c %s",
                              shQuote(
                                paste(
                                  "for f in",
                                  paste(shQuote(files), collapse = " "),
                                  "; do",
                                  "awk 'NR % 2 == 0' \"$f\" | sort -u >",
                                  "$f.kmers;",
                                  "done"
                                )
                              ))
              system(cmd1)

              # 2) Combine & filter shared kmers
              kmers_files <- paste0(files, ".kmers")

              # build the full pipeline in one string, then shQuote() it
              pipeline <- paste(
                "cat",
                paste(shQuote(kmers_files), collapse = " "),
                "|",
                "sort | uniq -c | awk '$1 >=",
                min_occurrence,
                "{print $2}' >",
                shQuote(shared_out)
              )

              cmd2 <- sprintf("bash -c %s", shQuote(pipeline))
              system(cmd2)

              # 3) Cleanup
              cmd3 <- sprintf("bash -c %s",
                              shQuote(paste("rm", paste(
                                shQuote(kmers_files), collapse = " "
                              ))))
              system(cmd3)

              message("  Wrote shared kmers: ", shared_out)
              out_files[grp] <- shared_out
            }

            out_files
          })
