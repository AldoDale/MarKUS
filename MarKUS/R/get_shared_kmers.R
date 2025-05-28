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

            # first, produce the shared_kmers per group as before
            shared_files <- setNames(character(length(grsplit)), names(grsplit))
            for (grp in names(grsplit)) {
              dfg        <- grsplit[[grp]]
              files      <- dfg$kmer_txt
              dir0       <- dirname(files[1])
              shared_out <- file.path(dir0, paste0(grp, "_shared_kmers.txt"))

              message("Working on group: ", grp)

              # 1) per-file unique and combine by min_occurrence
              cmd1 <- sprintf("bash -c %s",
                              shQuote(paste(
                                "for f in", paste(shQuote(files), collapse=" "), "; do",
                                "awk 'NR % 2 == 0' \"$f\" | sort -u > \"$f.kmers\";",
                                "done"
                              )))
              system(cmd1)

              kmers_files <- paste0(files, ".kmers")
              cmd2 <- paste(
                "cat", paste(shQuote(kmers_files), collapse = " "),
                "| sort | uniq -c | awk '$1 >=", min_occurrence, "{print $2}' >",
                shQuote(shared_out)
              )
              system(sprintf("bash -c %s", shQuote(cmd2)))

              # cleanup per-file .kmers immediately
              unlink(kmers_files)

              shared_files[grp] <- shared_out
            }

            # now do unique-to-group filtering in bash
            unique_files <- setNames(character(length(shared_files)), names(shared_files))
            for (grp in names(shared_files)) {
              this_f    <- shared_files[grp]
              others    <- shared_files[names(shared_files) != grp]
              unique_out <- sub("_shared_kmers\\.txt$", "_unique_kmers.txt", this_f)

              # build a grep pipeline: start from this_f, then subtract each other
              # -F: fixed, -x: whole-line, -v: invert, -f: patterns-from-file
              grep_cmds <- paste0(
                "grep -F -xv -f ", shQuote(others),
                collapse = " | "
              )
              full_cmd <- sprintf(
                "bash -c %s",
                shQuote(sprintf("cat %s | %s > %s",
                                this_f,
                                grep_cmds,
                                unique_out))
              )

              system(full_cmd)
              message("  Wrote UNIQUE kmers: ", unique_out)
              unique_files[grp] <- unique_out
            }

            unique_files
          })
