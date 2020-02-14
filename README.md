# kmer-tools

#### Introduction
The *canonical k-mer* will be the first alphabetically among the kmer and its reverse compliment, leaving a grand total of **32,768** possible 8-mers.

#### Kmer Lookup Table
A kmer lookup table will be created using an indexed list of alphabetically sorted kmers. 

#### Storage of kmer presence
Kmers can be represented as a relative frequency `N / (seqLength/k)` or by a boolean (present/absent) value.  To keep the storage space of each sequence in a kmer matrix low, only one byte of data will be used per 8-mer which results in a total of ~32KB per sequence.

##### ... as Relative Frequency
Relative frequency will be represented using *printable* ASCII values ranging from **33** to **126** ( ! to ~ ).
This gives the possible representation for frequencies between 0 and 94%. ...Assuming no 8-mer is in abundance greater than 94%.

|      | k0 | k1 | ... | k32767 |
|------|----|----|-----|--------|
| Seq1 | **!**  | c  | ... | 0      |
| Seq2 | a  | 9  | ... | x      |
| Seq3 | +  | x  | ... | **!**      |

##### ... as Boolean
If relative frequency data is available, boolean values will be determined based on NOT-PRESENT = ! and PRESENT is anything other than !.  Otherwise, the ASCII characters "1" *(49)* and "0" *(48)* will represent present and not-present respectively.

|      | k0 | k1 | ... | k32767 |
|------|----|----|-----|--------|
| Seq1 | **0**  | 1  | ... | 1      |
| Seq2 | 1  | 1  | ... | 1      |
| Seq3 | 1  | 1  | ... | **0**      |


