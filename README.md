# PSL GWAS
Eli Pandolfo

## Description
This project aims to perform a genome-wide association study (GWAS) on microbial
genomes (specifically E. coli), using Probabalistic Soft Logic (PSL).

A GWAS is a computational method of mapping specific attributes of a genotype
(such as SNPs, indels, or genes) to specific attributes of a phenotype (in this
case, resistance to different types of antibiotics). Formulated as a machine
learning problem, our model uses Bayesian inference to predict which genotype
attributes are associated with resistance to certain antibiotics, given
the presence or absence of that genotype attribute in a sample that either has
or lacks resistance to those antibiotics. We use kmers of length 30 to represent
genotype attributes, since they capture both point mutations and entire genes,
and are usable with unaligned genomes (we have only contigs of length 1000 to
100000 bps).

1. The first stage in the process is converting a directory full of FASTQ contig
files into a python nested dictionary, and serializing it with pickle. I
already did this. The associated program is pickle_input.py.

1. The second stage is generating a list of unique kmers in the dataset. I used
DSK for that.

1. The third stage is pruning that list down to kmers that appear at least
10 times, since we can't get any significant results on kmers that appear
rarely. The associated program is reduce_input.py.

1. The fourth stage is assembling a database that maps kmers to samples in which
those kmers appear. This allows us to have a contains(sample, kmer) lookup
table. The associated program is create_kmer_db.py. This is extremely
memory-intensive and you can't use a laptop to do it unless you have the power
of resurrection so you can see if it's finished after you die of boredom and
old age a few times. It scales non-linearly with available memory because of
the glory of O(1) hash table lookups.

1. The fifth stage is creating a reduced version of that database for PSL. We need
the full thing for the post-psl analysis at the end, but PSL does not need
to know where in the genome a kmer occurs. The associated program is
clean_kmer_db.py.

1. The sixth stage is converting the reduced database to the exact form PSL likes
as an input. The associated program is TBD.

1. The seventh stage is converting the table that logs which samples are resistant
to which antibiotics to a form that PSL likes as an input.

1. The eighth stage is running PSL.

1. The ninth stage is examining the subset of kmers that PSL flags,
and associating them with specific genes, indels, and SNPs.

1. The tenth stage is comparing those attributes to existing datasets,
to determine the accuracy of our method and to hopefully find new genes,
indels, and SNPs that contribute to resistance.


## Installation
You need python 3.6 or later, but there's no reason not to use
the latest version of python3. 

This project requires no external dependencies as of yet. I imagine
it will be easiest to use pandas for converting the tsv into a psl data
file form.

The only other thing to note here is that you should have a subdirectory
called `data/` in the root of this repo. That holds all the data files.

## Usage
Since the output of each stage is saved, you only need to do each stage once.
We are on stage 6 for our dataset. Since we haven't written the script
that does this, there are no usage instructions yet.

## Notes
General: identify the kmers that occur almost uniquely in
resistant samples and don't occur in non-resistant samples.
```
contains(s, k) & resistance(s, a) -> resistance(k, a)
contains(s, k) & !resistance(s, a) -> !resistance(k, a)
```

Extension: epistatic reactions between kmers.
```
contains(s, k1, k2) & resistance(s, a) & !resistance(k1, a) &
    !resistance(k2, a) -> positiveepistatic(k1, k2, a)
contains(s, k1, k2) & !resistance(s, a) & resistance(k1, a) &
    resistance(k2, a) -> negativeepistatic(k1, k2, a)
```
I'm not sure if PSL can handle this. Lotta ground rules. May be
ways to reduce tho.

More extension: increase to reactions between 3 kmers. Pipe dream.

The full schema for this repository, were you going to run the entire pipeline,
would look something like this. Note that I assume if you have a local DSK
binary, it will be somewhere else.

```
/........................................script associated with each stage of the pipeline
    README.md
    pickle_input.py......................stage 1
    reduce_input.py......................stage 2
    create_kmer_db.py....................stage 4
    clean_kmer_db.py.....................stage 5
    TDB..................................stage 6
    TBD..................................stage 7
    TBD..................................stage 8
    TBD..................................stage 9
    TBD..................................stage 10
    data/................................output of each stage of the pipeline
        raw.pickle.......................stage 1
        unique_kmers.txt.................stage 2
        unique_kmers_reduced.txt.........stage 3
        kmer_sample_map.txt..............stage 4
        kmer_sample_map_reduced.txt......stage 5
        TBD..............................stage 6
        TBD..............................stage 7
        TBD..............................stage 8
        TBD..............................stage 9
        TBD..............................stage 10
    contigs/
        sample1.fa
        sample2.fa
        ...
```

