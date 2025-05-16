# MarKUS: Marker K-mers for Untargeted Sequencing

MarKUS provides a fast and efficient workflow for the extraction, count, and filtering of k-mers from .fastq, .fasta, and .fastq.gz files.

## Installation

```r
remotes::install_github("AldoDale/MarKUS", subdir = "MarKUS")
```

## How to use


### Merge fastq files

If the sequencing data is split in first and second strand (or forward and reverse strand), the files should be first merged

```r

merged_files <- merge_fastq(path = ".", pattern = "*.fastq.1.gz")

#arguments
#  -path = path to the directory containing the files
#  -pattern = pattern to recognize the first (or forward) strand files

```


```r

count_kmers(x, kmer_size = 20, threads = 10)

```


```r


```
