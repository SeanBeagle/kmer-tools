# kmer-tools

#### Introduction
The *canonical k-mer* will be the first alphabetically among the kmer and its reverse compliment, resulting in a total of **32,768** possible 8-mers.

#### Kmer Lookup Table
A kmer lookup table will be created using an indexed list of alphabetically sorted kmers. 

|index| kmer|
|----|----|
|0|AAAAAAAA|
|1|AAAAAAAC|
|...|...|
|32767| NNNNNNNN |

#### Kmer Matrix
Kmers can be represented as a relative frequency `count / (seqLength-k)` or by a boolean `0 | 1` value. Each 8-mer will be allocated 1 Byte of data per sequence, resulting in a total of ~32KB space used per sequence.

##### ... as Relative Frequency
Relative frequency will be represented using the *94 printable* ASCII values ranging from **33** to **126** ( ! to ~ ).
This allows the possible values of relative frequencies to be between 0 and 94%. This assumes no 8-mer will be in abundance greater than 94%.

|      | k0 | k1 | ... | k32767 |
|------|----|----|-----|--------|
| Seq1 | **!**  | c  | ... | 0      |
| Seq2 | a  | 9  | ... | x      |
| Seq3 | +  | x  | ... | **!**      |

##### ... as Boolean
If relative frequency data is available, boolean states can be determined based on *ABSENT* as ` ! ` and *PRESENT* as anything other than ` ! `.  Otherwise, the ASCII characters "1" *(49)* and "0" *(48)* will represent PRESENT and ABSENT respectively.

|      | k0 | k1 | ... | k32767 |
|------|----|----|-----|--------|
| Seq1 | **0**  | 1  | ... | 1      |
| Seq2 | 1  | 1  | ... | 1      |
| Seq3 | 1  | 1  | ... | **0**      |

#### Machine Learning
... insert magic here!
