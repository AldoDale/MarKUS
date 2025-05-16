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
#  - path = path to the directory containing the files
#  - pattern = pattern to recognize the first (or forward) strand files

```

### Produce and count k-mers

To use this function, Jellyfish (https://github.com/gmarcais/Jellyfish) must be installed.

```r
kmers <- count_kmers(x, kmer_size = 20, threads = 10)

#arguments
#  - x = data.frame with a "sample" column with sample names and a "merged_path" column
         with the path to the fastq files. If merge_fastq() was used, the output of the
         function can be used as input of count_kmers()
#  - kmer_size = the desired size of k-mers
#  - threads = the number of threads to use

```

### Extract group-exclusive k-mers

This function is used to detect k-mers exclusive of one group of samples. A minimum prevalence value can be set.

```r
kmers <- get_shared_kmers(x, data, group = “origin”, min_occurrence = 4)


#arguments
#  - x = data.frame with a "sample" column with sample names and a "merged_path" column
         with the path to the fastq files. If merge_fastq() was used, the output of the
         function can be used as input of count_kmers()
#  - data = the desired size of k-mers
#  - group = the number of threads to use
#  - min_occurrence = the number of threads to use

```

```r


```
