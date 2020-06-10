Overview
#################

PSL-GWAS uses variable-length *k*-mers as genetic elements. It first enumerates all
*k*-mers present in all samples, and performs an initial filtering step to
remove *k*-mer, phenotype pairs not strongly correlated with each other (with a
customizable correlation threshold). It also computes
a similarity score for each pair of samples.

It then creates a Probabilistic Soft Logic (PSL) model, which outputs a score between
0 and 1 (inclusive) for each *k*-mer, phenotype pair. The PSL model implements a
kmowledge graph, where entities (samples, *k*-mers, and phenotypes) are nodes and
associations between those entities
(samples are similar, sample contains *k*-mer, sample displays phenotype,
phenotypes are similar, *k*-mer causes phenotype)
are edges. The PSL model efficiently performs convex optimization to find the
best satisfying assignment for the '*k*-mer causes phenotype' associations.

The PSL model is composed of first-order logical rules, which are customizable.
The base model, which generalizes across phenotypes and species, uses
similarity and dissimilarity between samples to model the bacterial population
structure without performing PCR, or MDS, and without making any assumptions
of independence among entities.

After running the association test, PSL-GWAS ranks the *k*-mer, phenotype
pairs by their confidence scores. The highest-scoring *k*-mers can be mapped
to an annotated reference genome and/or subjected to further analysis.

