#!/bin/bash
###############################################################################
# run.sh
# This script runs a GWAS. It checks for the existence of all preprocessed
# files used in the association test. If they do not exist, it runs the
# preprocessing pipeline to create them. It then invokes the PSL model to
# perform the association test. Finally, it invokes the postprocessing
# pipeline, and logs the settings used in the GWAS.
###############################################################################
ARGS="$@"
# define colors for pretty printing
RED="\033[0;31m"
GRN="\033[0;32m"
NC="\033[0m"
project=0
sample=0
pheno=0

check_project() {
    if [[ -d "$project" ]]; then
        echo "Found project $project"
    else
        echo "Could not find project $project."
        echo "Did you create it with the startproject script?"
        echo "Aborting."
        exit 1
    fi
}

check_sample() {
    sample_file="$project/data/raw/$sample"
    if [[ -e $sample_file ]]; then
        echo "Found sample file: $sample_file"
    else
        echo "Could not find sample file: $sample_file."
        echo "Are you sure this file exists?"
        echo "Aborting."
        exit 1
    fi
}

check_pheno() {
    pheno_file="$project/data/raw/$pheno"
    if [[ -e $pheno_file ]]; then
        echo "Found pheno file: $pheno_file"
    else
        echo "Could not find pheno file: $pheno_file."
        echo "Are you sure this file exists?"
        echo "Aborting."
        exit 1
    fi
}

check_usage() {
    if [[ "--" == $2 ]]; then
        echo "No $1 given"
        exit 1
    fi
}

set -- "$@" "--"

# extract options and their arguments into variables.
while (( "$#" )) ; do
    case "$1" in
        --project)
            check_usage "project name" $2
            project=$2
            shift 2;;
        --sample)
            check_usage "sample file name" $2
            sample=$2
            shift 2;;
        --pheno)
            check_usage "pheno file name" $2
            pheno=$2
            shift 2;;
        -d|--debug) shift;;
        --) shift;;
        *) echo "Invalid param" ; exit 1;;
    esac
done
check_project
check_sample
check_pheno

# define variables
OPATH="$project/data/postprocessed"
PPATH="$project/data/preprocessed"
RPATH="$project/data/raw"
Y="$GRN\xE2\x9C\x93$NC"
N="$RED\xE2\x9C\x97$NC"
postprocessed=0
preprocessed=0
raw=0

# check for postprocessed files
echo "Checking for data files..."
if [[ -r "$OPATH/scored_kmers.txt" ]] \
&& [[ -r "$OPATH/scored_kmers.fsa" ]]; then #TODO: add logs and meta check here
    postprocessed=1
    echo -e "\t$Y postprocessed"
else
    echo -e "\t$N postprocessed"
fi
# check for preprocessed files
if [[ -r "$PPATH/contains_sample_unitig.txt" ]] \
&& [[ -r "$PPATH/value_kmer_pheno.txt" ]] \
&& [[ -r "$PPATH/value_sample_pheno.txt" ]] \
&& [[ -r "$PPATH/similar_pheno_pheno.txt" ]] \
&& [[ -r "$PPATH/similar_sample_sample.txt" ]] \
&& [[ -r "$PPATH/pheno_int_map.pkl" ]] \
&& [[ -r "$PPATH/sample_int_map.pkl" ]] \
&& [[ -r "$PPATH/unitig_int_map.pkl" ]] \
&& [[ -r "$PPATH/unitig_sample_map.txt" ]] \
&& [[ -r "$PPATH/sample_similarities.tsv" ]] \
&& [[ -r "$PPATH/unique_kmers.txt" ]]; then
    preprocessed=1
    echo -e "\t$Y preprocessed"
else
    echo -e "\t$N preprocessed"
fi
# check for raw data
if [[ -r "$RPATH/$sample" ]] \
&& [[ -r "$RPATH/$pheno" ]]; then
    raw=1
    echo -e "\t$Y raw"
else
    echo -e "\t$N raw"
fi
echo

# runs preprocessing pipeline
run_preprocess() {
    echo "Running preprocessing pipeline"
    #./bin/preprocess.sh $ARGS
}

# runs psl assocations test
run_psl() {
    echo "Running association test"
    #./bin/run_psl.sh $ARGS
}

# runs postprocessing pipeline
run_postprocess() {
    echo "Running postprocessing pipeline"
    #./bin/postprocess.sh $ARGS
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

