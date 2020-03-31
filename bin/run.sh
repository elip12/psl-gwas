#!/bin/bash
###############################################################################
# run.sh
# This script runs a GWAS. It checks for the existence of all preprocessed
# files used in the association test. If they do not exist, it runs the
# preprocessing pipeline to create them. It then invokes the PSL model to
# perform the association test. Finally, it invokes the postprocessing
# pipeline, and logs the settings used in the GWAS.
###############################################################################
# define colors for pretty printing
RED='\033[0;31m'
GRN='\033[0;32m'
NC='\033[0m'

# ensure user inputs project name
if [[ "$#" -ne 1 ]]; then
    echo "Please specify a project name."
    exit 1
fi

# define variables
OPATH="$1/data/postprocessed"
PPATH="$1/data/preprocessed"
RPATH="$1/data/raw"
Y="$GRN\xE2\x9C\x93$NC"
N="$RED\xE2\x9C\x97$NC"
postprocessed=0
preprocessed=0
raw=0

# check for postprocessed files
echo "Checking for data files..."
if [[ -e "$OPATH/scored_kmers.txt" ]] \
&& [[ -e "$OPATH/scored_kmers.fsa" ]]; then #TODO: add logs and meta check here
    postprocessed=1
    echo -e "\t$Y postprocessed"
else
    echo -e "\t$N postprocessed"
fi
# check for preprocessed files
if [[ -e "$PPATH/contains_sample_kmer.txt" ]] \
&& [[ -e "$PPATH/resistance_kmer_class.txt" ]] \
&& [[ -e "$PPATH/resistance_sample_class.txt" ]] \
&& [[ -e "$PPATH/similar_pheno_pheno.txt" ]] \
&& [[ -e "$PPATH/similar_sample_sample.txt" ]] \
&& [[ -e "$PPATH/pheno_int_map.pkl" ]] \
&& [[ -e "$PPATH/sample_int_map.pkl" ]] \
&& [[ -e "$PPATH/unitig_int_map.pkl" ]] \
&& [[ -e "$PPATH/unitig_sample_map.txt" ]] \
&& [[ -e "$PPATH/similarities.tsv" ]] \
&& [[ -e "$PPATH/unique_kmers.tsv" ]]; then
    preprocessed=1
    echo -e "\t$Y preprocessed"
else
    echo -e "\t$N preprocessed"
fi
# check for raw data
if [[ -e "$RPATH/samples.tsv" ]] \
&& [[ -e "$RPATH/phenos.tsv" ]]; then
    raw=1
    echo -e "\t$Y raw"
else
    echo -e "\t$N raw"
fi
echo

# runs preprocessing pipeline
run_preprocess() {
    echo "Running preprocessing pipeline"
    ./bin/preprocess.sh $1
}

# runs psl assocations test
run_psl() {
    echo "Running association test"
    ./bin/prep_psl_data.sh $1
    ./bin/run_psl.sh $1
}

# runs postprocessing pipeline
run_postprocess() {
    echo "Running postprocessing pipeline"
    ./bin/postprocess.sh $1
}

# run required pipelines idempotently
if [[ $postprocessed -eq 1 ]]; then
    echo "Postprocessed data files found; nothing to do."
    echo "Exiting."
    exit 0
elif [[ $preprocessed -eq 1 ]]; then
    echo "Preprocessed data files found."
    run_psl && run_postprocess
    echo "Done."
    exit 0
elif [[ $raw -eq 1 ]]; then
    echo "Raw data files found."
    run_preprocess && run_psl && run_postprocess
    echo "Done."
    exit 0
else
    echo "Raw data files not found; nothing to do"
    echo "Exiting."
    exit 0
fi

