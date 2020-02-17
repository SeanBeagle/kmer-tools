# kmer-counter-master & kmerizer_to_csv

#### Introduction
The *canonical k-mer* will be the first alphabetically among the kmer and its reverse compliment, resulting in a total of **32,768** possible 8-mers.

#### Install Kmer-Counter-Master
Unzip kmer-counter-master prior to identifying kmer frequencies. 

``` 
cd scripts/
unzip kmer-counter-master.zip
cd kmer-counter-master/
make
cd ../..
```

#### Running Kmerizer_to_csv
Kmerizer_to_csv is a python 3 argparse script with required input of a fasta file (no checks currently in place). The default kmer size is 8 (hard coded with columns to keep, can change) and a default output is provided. The .csv extension is automatically applied. 

```
scripts/kmerizer_to_csv.py -f [input_file] -k 8 -o [output_file_name]

```