Overview
#################

PSL-GWAS uses unitigs as genetic elements. It first enumerates all
unitigs present in all samples, and performs an initial filtering step to
remove all unitig-pheno pairs not associated with each other. It also generates
a similarity value for each pair of samples.

It then creates a Probabilistic Soft Logic (PSL) model, which outputs a score between
0 and 1 (inclusive) for each unitig-phenotype pair. The PSL model implements a
kmowledge graph, where entities (samples, unitigs, and phenotypes) are nodes and
associations between those entities
(samples are similar, sample contains unitig, sample displays phenotype,
phenotypes are similar, unitig causes phenotype)
are edges. The PSL model efficiently performs convex optimization to find the
best satisfying assignment for the 'unitig causes phenotype' associations.

The PSL model is composed of first-order logical rules, which are customizable.
The base model, which generalizes across phenotypes and species, uses
similarity and dissimilarity between samples to model the bacterial population
structure without performing PCR, or MDS, and without making any assumptions
of independence among entities.

After running the association test, PSL-GWAS extracts the unitig-phenotype
pairs with confidence above a customizable threshold, to be mapped to a
reference genome and/or subjected to further analysis.

