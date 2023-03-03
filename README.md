# *aardvark*

## Dependencies

* A c++17 supporting compiler

## Installation

```bash
git clone https://github.com/holstegelab/aardvark.git && cd aardvark
mkdir build
make
```

## Cyclical-shift-aware kmer-counting

```bash
//count kmers of length 21
aardvark kc -k 21 *.fa > counts.tsv

//count all kmers from length 3 to 32
aardvark kc -k 3 -m 32 *.fa > counts.tsv

```
Note: maximum kmer-length is 32.

## Maximum-sequence masking

```bash
//mask all subsequences in the given FASTA files that are present in 'seqs_to_mask.fa'
aardvark mk -m seqs_to_mask.fa *.fa > masked.fa

```