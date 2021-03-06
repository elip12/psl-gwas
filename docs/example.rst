PSL-GWAS Example
################

We provide example data to allow you to replicate the results of `Pandolfo et al. 2020
<nothing yet>`_.

#. Ensure postgres is installed and running, and create a new database your default user can access.
#. Clone the repo and cd into it.
#. Download the example data using the ``fetch_example`` script.
#. Run the GWAS.

.. code-block:: bash

    git clone https://github.com/elip12/psl-gwas.git
    cd psl-gwas
    ./bin/fetch_example.sh
    ./bin/run.sh --project example --sample samples.tsv --pheno phenos.tsv --truth truths.fa -p --postgres <DATABASE>

In ``example/data/postprocessed/`` you will see one file for each phenotype, holding the highest-ranking *k*-mers.

To evaluate the results, you can use the included evaluation script. The first argument
is the path to the truth data, the second is the path to output *k*-mers, and the third
is the number of *k*-mers to use when evaluating.

.. code-block:: bash

    python3 psl-gwas/extra/evaluate.py example/data/raw/truths.fa example/data/postprocessed/psl.tobramycin 10

This invocation will compute and print the precision, hits, and recirpocal rank of the first hit, as well
as printing out the annotations of the specific truth sequences identified.
