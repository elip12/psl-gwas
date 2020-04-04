#!/bin/bash
###############################################################################
##  run.sh
##  This script runs a GWAS. It checks for the existence of all preprocessed
##  files used in the association test. If they do not exist, it runs the
##  preprocessing pipeline to create them. It then invokes the PSL model to
##  perform the association test. Finally, it invokes the postprocessing
##  pipeline.
###############################################################################
# define colors for pretty printing
RED="\033[0;31m"
GRN="\033[0;32m"
NC="\033[0m"
project=0
sample=0
pheno=0

# checks project directory exists
check_project() {
    if [[ -d "$project" ]]; then
        echo "Found project: $project"
    else
        echo "Could not find project $project"
        echo "Did you create it with the startproject script?"
        echo "Aborting."
        exit 1
    fi
}

# checks user-given sample file exists
check_sample() {
    sample_file="$project/data/raw/$sample"
    if [[ -e $sample_file ]]; then
        echo "Found sample file: $sample_file"
    else
        echo "Could not find sample file: $sample_file"
        echo "Are you sure this file exists?"
        echo "Aborting."
        exit 1
    fi
}

# checks user-given pheno file exists
check_pheno() {
    pheno_file="$project/data/raw/$pheno"
    if [[ -e $pheno_file ]]; then
        echo "Found pheno file: $pheno_file"
    else
        echo "Could not find pheno file: $pheno_file"
        echo "Are you sure this file exists?"
        echo "Aborting."
        exit 1
    fi
}

# ensures users give param value when one is required
check_usage() {
    if [[ "--" == $2 ]]; then
        echo "$1 requires an argument"
        exit 1
    fi
}

# runs preprocessing pipeline
run_preprocess() {
    echo "Running preprocessing pipeline"
    ./bin/preprocess.sh ${pre_opts[@]}
}

# runs psl assocation test
run_psl() {
    echo "Running association test"
    ./bin/run_psl.sh ${psl_opts[@]}
}

# runs postprocessing pipeline
run_postprocess() {
    echo "Running postprocessing pipeline"
    ./bin/postprocess.sh ${pre_opts[@]}
}

trap exit SIGINT 
set -- "$@" "--"
psl_opts=()
pre_opts=()
# parses command line options
while (( "$#" )) ; do
    case "$1" in
        --project)
            check_usage "--project" $2
            project=$2
            psl_opts+=("--project" "$2")
            pre_opts+=("--project" "$2")
            shift 2;;
        --sample)
            check_usage "--sample" $2
            sample=$2
            pre_opts+=("--sample" "$2")
            shift 2;;
        --pheno)
            check_usage "--pheno" $2
            pheno=$2
            pre_opts+=("--pheno" "$2")
            shift 2;;
        -d|--debug)
            pre_opts+=("--debug")
            shift;;
        --threads)
            pre_opts+=("--threads" "$2")
            check_usage "--threads" $2
            shift 2;;
        --mem)
            check_usage "--mem" $2
            psl_opts+=("--mem" "$2")
            pre_opts+=("--mem" "$2")
            shift 2;;
        -k|--k)
            check_usage "-k/--k" $2
            pre_opts+=("--k" "$2")
            shift 2;;
        --upperfreq)
            check_usage "--upperfreq" $2 
            pre_opts+=("--upperfreq" "$2")
            shift 2;;
        --lowerfreq)
            check_usage "--lowerfreq" $2
            pre_opts+=("--lowerfreq" "$2")
            shift 2;;
        --thresh)
            check_usage "--thresh" $2 
            pre_opts+=("--thresh" "$2")
            shift 2;;
        --postgres)
            check_usage "--postgres" $2
            psl_opts+=("${psl_opts[@]}" "--postgres" "$2")
            shift 2;;
        -p|--param)
            pre_opts+=("--param")
            shift;;
        --) shift;;
        *) echo "Invalid argument: $1" ; exit 1;;
    esac
done
check_project

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
if [[ -r "$OPATH/scored_unitigs.txt" ]] \
&& [[ -r "$OPATH/scored_unitigs.fsa" ]] \
&& [[ -r "$OPATH/UNITIGPHENO.txt" ]]; then
    postprocessed=1
    echo -e "$Y postprocessed"
else
    echo -e "$N postprocessed"
fi
# check for preprocessed files
if [[ -r "$PPATH/contains_sample_unitig.txt" ]] \
&& [[ -r "$PPATH/value_unitig_pheno.txt" ]] \
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
    echo -e "$Y preprocessed"
else
    echo -e "$N preprocessed"
fi
# check for raw data
if [[ -r "$RPATH/$sample" ]] \
&& [[ -r "$RPATH/$pheno" ]]; then
    raw=1
    echo -e "$Y raw"
else
    echo -e "$N raw"
fi

# run required pipelines idempotently
if [[ $postprocessed -eq 1 ]]; then
    echo "Postprocessed data files found; nothing to do."
    echo "Exiting."
    exit 0
elif [[ $preprocessed -eq 1 ]]; then
    if ! [[ -r "$OPATH/UNITIGPHENO.txt" ]]; then
        run_psl
    fi
    run_postprocess
    echo "Done."
    exit 0
elif [[ $raw -eq 1 ]]; then
    run_preprocess && run_psl && run_postprocess
    echo "Done."
    exit 0
else
    echo "Raw data files not found; nothing to do."
    echo "Exiting."
    exit 0
fi
