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
    if [[ -e $sample_file ]];
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
    if [[ -e $pheno_file ]];
        echo "Found pheno file: $pheno_file"
    else
        echo "Could not find pheno file: $pheno_file."
        echo "Are you sure this file exists?"
        echo "Aborting."
        exit 1
    fi
}

# ensure user inputs necessary args
TEMP=`getopt -o d::k::\
--long project:,samples:,phenos:,threads::,mem::,upperfreq::,lowerfreq::,\
thresh::,param::,debug:: -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
    case "$1" in
        --project)
            project=$2
            check_project ; shift 2;;
        --sample)
            sample=$2 ; shift 2;;
        --pheno)
            pheno=$2 ; shift 2;;
        -d|--debug) ; shift;;
        --) shift 2; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

check_sample
check_pheno

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
if [[ -e "$PPATH/contains_sample_unitig.txt" ]] \
&& [[ -e "$PPATH/value_kmer_pheno.txt" ]] \
&& [[ -e "$PPATH/value_sample_pheno.txt" ]] \
&& [[ -e "$PPATH/similar_pheno_pheno.txt" ]] \
&& [[ -e "$PPATH/similar_sample_sample.txt" ]] \
&& [[ -e "$PPATH/pheno_int_map.pkl" ]] \
&& [[ -e "$PPATH/sample_int_map.pkl" ]] \
&& [[ -e "$PPATH/unitig_int_map.pkl" ]] \
&& [[ -e "$PPATH/unitig_sample_map.txt" ]] \
&& [[ -e "$PPATH/sample_similarities.tsv" ]] \
&& [[ -e "$PPATH/unique_kmers.txt" ]]; then
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
    ./bin/preprocess.sh $ARGS
}

# runs psl assocations test
run_psl() {
    echo "Running association test"
    ./bin/prep_psl_data.sh $ARGS
    ./bin/run_psl.sh $ARGS
}

# runs postprocessing pipeline
run_postprocess() {
    echo "Running postprocessing pipeline"
    ./bin/postprocess.sh $ARGS
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

