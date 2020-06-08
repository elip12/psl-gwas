Installation
############

Required dependencies:

* Bash
* Python3
* Java

Recommended dependencies:

* Postgres

We recommend using a python3 virtual environment (or some equivalent like
`pyenv`_).

Installation:

.. code-block:: bash

    git clone https://github.com/elip12/psl-gwas.git
    cd psl-gwas
    ./bin/startproject.sh <project>

This will create a project directory with the given project name,
copy all necessary files into it, and install all necessary python packages.

The project directory has the following structure::

    <project>/
    ....gwas.psl            defines PSL model
    ....gwas.data           defines PSL data files
    ....parameters.yaml     optional, replaces some command-line params
    ....data/               holds data files
    ........raw/            holds raw data files
    ........preprocessed/   holds preprocessed data files
    ........postprocessed/  holds postprocessed data files

.. _pyenv: https://github.com/pyenv/pyenv

