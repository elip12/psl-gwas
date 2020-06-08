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

