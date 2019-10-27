What we have:
    - a file containing the ~18 million unique kmers in the dataset that occur
    10 or more times
    - the unaligned genome of each of the ~355 samples
    - a mapping between samples and resistance
What we need:
    - a mapping between all samples and all unique kmers
Why we can't get this:
    - it takes a fucking long time: 18M * 355 * 3.5M comparisons

## Notes

input.txt is a file containing each unique kmer of size 30 in the dataset,
which i got by running dsk.

raw.pickle is the entire input dataset as a dict of dicts of lists of strings.

The idea is to map each kmer to the sample it comes from, then use the
knowledge of which samples are resistant to infer which kmers cause resistance.

contains(s, k) & resistance(s) -> resistance(k)
contains(s, k) & !resistance(s) -> !resistance(k)

Generalization: identify the kmers that occur almost uniquely in
resistant samples and don't occur in non-resistant samples.

Extension: epistatic reactions between kmers



