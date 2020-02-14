# kmer-tools

#### Introduction
The *canonical k-mer* will be the first alphabetically of the kmer and its reverse compliment.
The total number of 8-mers is equal to `(4^8)/2` = **32,768**.

#### Kmer Lookup Table
A kmer lookup table will be sorted alphabetically and each kmer will be assigned an index (0-32767). 

#### Storage of kmer presence
Kmers can be identified by their frequency relative to the total possible 8-mers `len(seq)/8` or by a boolean (present/absent) value.  To keep the footprint of each sequence object in a kmer matrix low, one byte of data will be used per 8-mer which results in a total of ~32KB per sequence.

##### ... as Relative Frequency
Relative frequency will be represented using *printable* ascii values ranging from **33** to **126** ( ! to ~ ).
This gives the possible representation for frequencies between 0 and 94%. ...Assuming no 8-mer is in abundance greater than 94%.
##### ... as Boolean
If relative frequency data is available, boolean values will be determined based on NOT-PRESENT = ! and PRESENT is anything other than !.


|-|-|-|-|
| | K1 | K2| K3 |
OBJ1 | ! | a | ~ |
OBJ2 | ! | x | 9 | 
OBJ3 | ! | a | i |


#TODO: 
1. Sean make product class
