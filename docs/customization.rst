Customizing Your Model
######################

The base PSL model (defined in ``<project>/gwas.psl``) has the following rules:

1. Contains(S, K) & SamplePheno(S, P) >> KmerPheno(K, P)
    If a sample contains a *k*-mer and displays a pheno, there is evidence that
    *k*-mer causes that pheno.
2. Contains(S, K) & !SamplePheno(S, P) >> !KmerPheno(K, P)
    If a sample contains a *k*-mer and does not display a pheno, there is
    evidence that *k*-mer does not cause that pheno.
3. SimilarPheno(P1, P2) & KmerPheno(K, P1) & (P1 != P2) >> KmerPheno(K, P2)
    If two phenos are similar, and a *k*-mer causes one of them, there is evidence
    it causes the other. Note that PSL-GWAS by default only considers two
    phenos 'similar' if they are very strongly correlated.
4. !SimilarPheno(P1, P2) & KmerPheno(K, P1) & (P1 != P2) >> !KmerPheno(K, P2) ^2
    Two dissimilar phenos are unlikely to be caused by the same kmer.
5. SimilarSample(S1, S2) & SamplePheno(S1, P) & !SamplePheno(S2, P) & Contains(S1, K) & !Contains(S2, K) & (S1 != S2) >> KmerPheno(K, P)
    If two samples are similar, and one contains a *k*-mer and displays a pheno and
    the other does neither, there is evidence that *k*-mer causes that pheno.
    This is more tailored toward finding causal SNPs. An extreme example of this could
    be two sister bacteria, one of which developed a mutation during mitosis
    that confers a phenotype.
6. DissimilarSample(S1, S2) & SamplePheno(S1, P) & SamplePheno(S2, P) & Contains(S1, K) & Contains(S2, K) & (S1 != S2) >> KmerPheno(K, P)
    If two samples are dissimilar, and both display a pheno and both contain the
    same *k*-mer, there is evidence that *k*-mer causes the pheno.
    This is more tailored toward finding causal genes. An extreme example of this could
    be two bacteria of two different species that share almost no chromosomal DNA but
    share a phenotype-causing plasmid.
7. !DissimilarSample(S1, S2) & SamplePheno(S1, P) & SamplePheno(S2, P) & Contains(S1, K) & Contains(S2, K) & (S1 != S2) >> !KmerPheno(K, P) ^2
    This is the opposite of the above rule; it serves to lower the confidence scores
    of kmers that appear in highly related samples. 
8. !KmerPheno(K, P)
    Most *k*-mers do not cause phenos.

You can add, delete, or reweight any rules you like. For example, if you suspect phenotypic
variation in your data is primarily driven by SNPs, you could increase the weight
of rule 5. If you know two phenotypes usually do not occur in the same sample,
you could add a rule formalizing that logic. Further, you can replace
the data files PSL-GWAS generates with your own; you could use something like
`ClonalFrameML`_ to generate a recombination-aware model of the population structure,
and use that as the data for the SimilarSample and DissimilarSample rules.

For more information, see `PSL rule specification`_.

.. _PSL rule specification: https://psl.linqs.org/wiki/master/Rule-Specification.html
.. _ClonalFrameML: https://github.com/xavierdidelot/clonalframeml
