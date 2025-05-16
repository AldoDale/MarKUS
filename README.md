# MarKUS: Marker K-mers for Untargeted Sequencing

MarKUS provides a fast and efficient workflow for the extraction, count, and filtering of k-mers from .fastq, .fasta, and .fastq.gz files.

## Installation

```r
# or from GitHub
remotes::install_github("AldoDale/MarKUS", subdir = "MarKUS")
```

## How to use

### If the sequencing data is split in first and second strand (or forward and reverse strand), the files should be first merged

```r

merged_files <- merge_fastq(path = ".", pattern = "*.fastq.1.gz")

```
###arguments
  -dsj
  -dfsfd

```r

merged_files <- merge_fastq(path = ".", pattern = "*.fastq.1.gz")

```


```r


```
