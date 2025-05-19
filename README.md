# MarKUS: Marker K-mers for Untargeted Sequencing


MarKUS provides a fast and efficient workflow for the extraction, count, and filtering of k-mers from .fastq, .fasta, and .fastq.gz files.

## Installation

```r
remotes::install_github("AldoDale/MarKUS", subdir = "MarKUS")
```
---

## How to use

An example dataset is available at https://github.com/AldoDale/MarKUS/tree/main/MarKUS/man/example_dataset.


### Merge fastq files

If the sequencing data is split in first and second strand (or forward and reverse strand), the files should be first merged

```r

merge_fastq(path, pattern)

#arguments

#  - path = path to the directory containing the files.

#  - pattern = pattern to recognize the first (or forward) strand files.
```

#### Returns: a data.frame
```r

merged_files <- merge_fastq(path = ".", pattern = "*.fastq.1.gz")

merged_files

#>       sample                     merged_path
#>1     Malta_1    ./path/to/file/Malta_1.fastq
#>2     Malta_2    ./path/to/file/Malta_2.fastq
#>3     Malta_3    ./path/to/file/Malta_3.fastq

```
---

### Produce and count k-mers

To use this function, Jellyfish (https://github.com/gmarcais/Jellyfish) must be installed.

```r
count_kmers(x, kmer_size, threads)


#arguments

#  - x = data.frame with a "sample" column with sample names and a "merged_path" column
#        with the path to the fastq files. If merge_fastq() was used, the output of the
#        function can be used as input of count_kmers().
         
#  - kmer_size = the desired size of k-mers.

#  - threads = the number of threads to use.

```
#### Returns: a data.frame
```r

kmers <- count_kmers(merged_files, kmer_size = 20, threads = 10)

kmers

#>       sample                             counts_db                                    kmer_txt
#>1     Malta_1    ./path/to/file/Malta_1_counts_R.jf    ./path/to/file/Malta_1_kmer_counts_R.txt
#>2     Malta_2    ./path/to/file/Malta_2_counts_R.jf    ./path/to/file/Malta_2_kmer_counts_R.txt
#>3     Malta_3    ./path/to/file/Malta_3_counts_R.jf    ./path/to/file/Malta_3_kmer_counts_R.txt

```
---

### Extract group-exclusive k-mers

This function is used to detect k-mers exclusive of one group of samples. A minimum prevalence value can be set.

```r
get_shared_kmers(x, data, group, min_occurrence)


#arguments

#  - x = data.frame with a "sample" column with sample names and a "kmer_txt" column
#        with the path to the k-mers .txt files. The dataframe is created by the function count_kmers().
         
#  - data = a data.frame with the metadata. It must have a "sample" column (identical to the x data.frame column).

#  - group = the name of the metadata column to use as grouping factor.

#  - min_occurrence = the minimum number of samples the k-mers must be present in to be kept.

```

#### Returns: a named vector
```r

shared_kmers <- get_shared_kmers(kmers, data, group = “origin”, min_occurrence = 4)


shared_kmers

#>                                       Malta                                      Portugal                                      Spain
#>   "./path/to/file/Malta_shared_kmers_R.txt" "../path/to/file/Portugal_shared_kmers_R.txt" "../path/to/file/Spain_shared_kmers_R.txt" 

```

---

###Calculate information content of k-mers

This function calculates, based on the frequency of nitrogenous bases, the information content of k-mers. To find the most informative markers, the user might want to filter the k-mers based on this index.

```r
calculate_shannon(x)


#arguments

# - x = a named vector. Names must indicate the group, values must be the path to the "<group>_shared_kmers_R.txt" file.

```

#### Returns: a list of lists of named vectors and plots.
```r
shannon <- calculate_shannon(shared_kmers)

shannon$values

#>$Malta
#>AAAAGTGGAA AAAGTAGTGT AAATAACAGC AAATTCGCCT AACCCATCTA AAGAAACAGC AAGAATATGA AAGACTTTGT AAGCCCTATC AATACACGAG AATCCAACGG AATCTACCGC 
#> 1.2954618  1.5709506  1.5709506  1.8954618  1.5219281  1.3709506  1.3709506  1.8464393  1.8464393  1.7609640  1.8464393  1.8464393  
#> 
#>$Portugal
#>AAAAACTGCA AAAATTTTGG AAACGCTAGT AAACTTGACT AAAGGACTCG AAATGCAGTG AAATGTACGT AAATTACGCT AACACAATCG AACACTCAAG AACCGACGAC AACGGATTGT 
#>  1.570951   1.521928   1.921928   1.846439   1.846439   1.846439   1.846439   1.846439   1.685475   1.685475   1.521928   1.895462
#>  
#>$Spain
#>AAAACATAAT AAACGCAAGA AAATTACATA AACCGATCGG AACTCTCCAG AAGAAATAAA AAGCAGACCG AAGTCGGACT AATCGTAGCT AATTAATGCG AATTTGTTCT ACACACTGCC
#> 1.1567796  1.3709506  1.2954618  1.8954618  1.8464393  0.9219281  1.5709506  1.9709506  1.9709506  1.8464393  1.5709506  1.6854753
 
shannon$plots

#>$Malta

#>$Portugal

#>$Spain
 

```
<p align="center">
  <img src="MarKUS/man/example_figures/malta.png" width="400" />
  <img src="MarKUS/man/example_figures/portugal.png" width="400" />
  <img src="MarKUS/man/example_figures/spain.png" width="400" />
</p>

### Filter k-mers based on information content


```r

filter_shannon_values(x,  threshold)

#arguments

# - x = a list of named vectors (e.g., shannon$values).
# - threshold = minimum shannon value of sequences to be retained.
```

##### Returns: a list of named vectors

```r

filtered_shannon <- filter_shannon_values(shannon$values,  threshold = 1.8)
filtered_shannon 

#>$Malta
#>AGTTACGCTG ATCCGATTGA ATTTACCGGC CAACGTTTGG CACGTGTTGA CTGGATCAGC GCTTACTCGA GGCTTCACTA GTCGTTCAAC GTGCACACTA TCTACTGGAA 
#>  1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951 
#>
#>$Portugal
#>ACACTGTAGC ACACTGTCGT ACCTTCAAGG AGGTACCCTG ATATGGCGCA ATCTGTGGAC ATGGACACTC ATGTCAGGAC ATTCAGAGTC CAATCGGCTG CAGTCATGCA 
#>  1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   
#>CGTACGTGCA GACTGTAACC GCTATCAGCA GGAAGTCTCC GGACTGCTAC TCTATGGACA TCTGGAACCA 
#>  1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951 
#>
#>$Spain
#>AAGTCGGACT AATCGTAGCT ACCAGTTTGA ACGCGACTAT AGAGCCTTTA AGCAGTTACG AGCTTATAGC ATCCTCTGGA CAAGGTTGCA CAGTGACTTC CCTGCTGTAA 
#>  1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   1.970951   
#>GAGTCATCAC GCCAGGTTAC GGAGTCATCA 
#>  1.970951   1.970951   1.970951 

```

### Filter k-mers based on the edit distance

The sequences which were found to be unique might be biased by mutations and sequencing errors (substitution, addition, deletion, or frameshift). The function filter_ed() allows to detect possible biases by calculating the number of edits that have to be done on a string (in this case a sequence) to make it identical to another string, and filter them based on a threshold. The function is based on the package stringdist (https://github.com/markvanderloo/stringdist) and uses its methods.

```r

filter_ed(x, threshold, method, chunk_size, PPARAM)

#arguments

# - x = a list of named vectors (e.g., shannon$values).
# - threshold = minimum number of edits for a sequence to be retained as unique. 
# - method = one of the methods inherited by stringdist package.
# - chunk_size = number of lines to be processed together. Useful when dealing with big datasets.
# - PPARAM = parameter for the parallelization of processes.
```

##### Returns: a list of named vectors

```r

filt_edit_dist <- filter_ed(filt_shann, threshold = 3, method = "lv", chunk_size = 1, BPPARAM = BiocParallel::MulticoreParam())

filt_edit_dist

#>$Malta
#>[1] "AGTTACGCTG" "ATCCGATTGA" "ATTTACCGGC" "CTGGATCAGC" "GCTTACTCGA" "GTCGTTCAAC" "TCTACTGGAA"
#>
#>$Portugal
#> [1] "ACCTTCAAGG" "AGGTACCCTG" "ATATGGCGCA" "ATGGACACTC" "CAATCGGCTG" "CAGTCGGTAC" "GACTGTAACC" #>"GCTATCAGCA"
#> [9] "GGAAGTCTCC" "GGACTGCTAC" "TCTGGAACCA"
#>
#>$Spain
#> [1] "ACCAGTTTGA" "ACGCGACTAT" "AGAGCCTTTA" "AGCAGTTACG" "AGCTTATAGC" "ATCCTCTGGA" "CAAGGTTGCA" #>"CAGTGACTTC"
#> [9] "CCTGCTGTAA" "CGCCGTATTA" "GCCAGGTTAC"
```

### Filtering based on the presence of repeated patterns

Sequences that passed the previous filtering steps might present repeated patterns which lower their informativeness in some context.

```r

filter_repeated_seqs(x, mode, pattern, min_repeats, kmer_length) 


#arguments

# - x = a list of named vectors (e.g., shannon$values).
# - mode = one of "pattern" or "any". With "pattern" the user needs to input the pattern to search for. "any" is used to search any             repeated sequence.
# - pattern = the pattern to search if mode = "pattern".
# - min_repeats = number of times a sequence or pattern needs to be present to exclude the k-mer.
# - kmer_length = the length of the sequence to search for if it is repeated if mode = "any".
```

##### Returns: a list of data.frames

```r

# By using 'mode = "pattern"' we can remove the k-mers which present the <pattern> <min_repeats> times.
# In this example, we want to remove sequences where "AT" is repeated at least 2 times.

filt_pattern <- <- filter_repeated_seqs(fe, mode = "pattern", pattern = "AT", min_repeats = 2, kmer_length = NULL) 

filt_pattern

#>$Malta
#>DataFrame with 7 rows and 2 columns
#>     sequence       group
#>  <character> <character>
#>1  AGTTACGCTG       Malta
#>2  ATCCGATTGA       Malta
#>3  ATTTACCGGC       Malta

#>$Portugal
#>DataFrame with 10 rows and 2 columns
#>      sequence       group
#>   <character> <character>
#>1   ACCTTCAAGG    Portugal
#>2   AGGTACCCTG    Portugal
#>3   ATGGACACTC    Portugal

#>$Spain
#>DataFrame with 11 rows and 2 columns
#>      sequence       group
#>   <character> <character>
#>1   ACCAGTTTGA       Spain
#>2   ACGCGACTAT       Spain
#>3   AGAGCCTTTA       Spain


#In this example we want to remove sequences with any k-mer of length 2 repeated at least 2 times.

filt_any <- filter_repeated_seqs(fe, mode = "any", pattern = NULL, min_repeats = 2, kmer_length = 2)

filt_any

#>$Malta
#>DataFrame with 3 rows and 2 columns
#>     sequence       group
#>  <character> <character>
#>1  AGTTACGCTG       Malta
#>2  ATTTACCGGC       Malta
#>3  CTGGATCAGC       Malta
#>
#>$Portugal
#>DataFrame with 4 rows and 2 columns
#>     sequence       group
#>  <character> <character>
#>1  ACCTTCAAGG    Portugal
#>2  AGGTACCCTG    Portugal
#>3  CAATCGGCTG    Portugal
#>4  TCTGGAACCA    Portugal
#>
#>$Spain
#>DataFrame with 3 rows and 2 columns
#>     sequence       group
#>  <character> <character>
#>1  ACCAGTTTGA       Spain
#>2  CAGTGACTTC       Spain
#>3  GCCAGGTTAC       Spain
```




