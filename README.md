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
#### Returns: a data.frame
```r


merged_files

#>       sample                     merged_path
#>1     Malta_1    ./path/to/file/Malta_1.fastq
#>2     Malta_2    ./path/to/file/Malta_2.fastq
#>3     Malta_3    ./path/to/file/Malta_3.fastq

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

###Extract group-exclusive k-mers

This function is used to detect k-mers exclusive of one group of samples. A minimum prevalence value can be set.

```r
shared_kmers <- get_shared_kmers(x, data, group = “origin”, min_occurrence = 4)


#arguments

#  - x = data.frame with a "sample" column with sample names and a "kmer_txt" column
         with the path to the k-mers .txt files. The dataframe is created by the function           count_kmers()
         
#  - data = a data.frame with the metadata. It must have a "sample" column (identical to            the x data.frame column).

#  - group = the name of the metadata column to use as grouping factor

#  - min_occurrence = the minimum number of samples the k-mers must be present in to be         kept

```

###Calculate information content of k-mers

This function calculates, based on the frequency of nitrogenous bases, the information content of k-mers. To find the most informative markers, the user might want to filter the k-mers based on this index.

```r
shannon <- calculate_shannon(x)


#arguments

# - x = a named vector. Names must indicate the group, values must be the path to the "<group>_shared_kmers_R.txt" file.

#output


```
