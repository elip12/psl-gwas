# PSL GWAS
Eli Pandolfo

## Description
This project aims to perform a genome-wide association study (GWAS) on microbial
genomes (specifically E. coli), using Probabalistic Soft Logic (PSL).

A GWAS is a computational method of mapping specific attributes of a genotype
(such as SNPs, indels, or genes) to specific attributes of a phenotype (in this
case, resistance to different types of antibiotics).

The first stage in the process is converting a directory full of FASTQ contig
files into a python nested dictionary, and serializing it with pickle. I
already did this. The associated program is pickle_input.py.

The second stage is generating a list of unique kmers in the dataset. I used
DSK for that.

The third stage is pruning that list down to kmers that appear at least
10 times, since we can't get any significant results on kmers that appear
rarely. The associated program is reduce_input.py.

The fourth stage is assembling a database that maps kmers to samples in which
those kmers appear. This allows us to have a contains(sample, kmer) lookup
table. The associated program is create_kmer_db.py. This is extremely
memory-intensive and you can't use a laptop to do it unless you have the power
of resurrection so you can see if it's finished after you die of boredom and
old age a few times. It scales non-linearly with available memory because of
the glory of O(1) hash table lookups.

The fifth stage is creating a reduced version of that database for PSL. We need
the full thing for the post-psl analysis at the end, but PSL does not need
to know where in the genome a kmer occurs. The associated program is
clean_kmer_db.py.

The sixth stage is converting the reduced database to the exact form PSL likes
as an input. The associated program is TBD.

The seventh stage is converting the table that logs which samples are resistant
to which antibiotics to a form that PSL likes as an input.

The eighth stage is running PSL.

The ninth stage is examining the subset of kmers that PSL flags,
and associating them with specific genes, indels, and SNPs.

The tenth stage is comparing those attributes to existing datasets,
to determine the accuracy of our method and to hopefully find new genes,
indels, and SNPs that contribute to resistance.


## Installation

## Usage

## Notes
General: identify the kmers that occur almost uniquely in
resistant samples and don't occur in non-resistant samples.
```
contains(s, k) & resistance(s) -> resistance(k)
contains(s, k) & !resistance(s) -> !resistance(k)
```

Extension: epistatic reactions between kmers.
```
contains(s, k1, k2) & resistance(s) & !resistance(k1) &
    !resistance(k2) -> positiveepistatic(k1, k2)
contains(s, k1, k2) & !resistance(s) & resistance(k1) &
    resistance(k2) -> negativeepistatic(k1, k2)
```
I'm not sure if PSL can handle this. Lotta ground rules. May be
ways to reduce tho.

More extension: increase to reactions between 3 kmers. Pipe dream.

