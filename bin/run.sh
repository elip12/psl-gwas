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

get_mem() {
    if [[ "$OSTYPE" == "linux-gnu" ]]; then
        mem=$(( $(free -g | awk '/^Mem:/{print $2}') * 19 / 20))
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        mem=$(( $(sysctl -n hw.memsize) * 19 / 20 / 2**30))
    fi
    echo $mem
}

get_cpus() {
    echo $(getconf _NPROCESSORS_ONLN)
}

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
proj=()
# parses command line options
while (( "$#" )) ; do
    case "$1" in
        --project)
            check_usage "--project" $2
            project=$2
            psl_opts+=("--project" "$2")
            proj+=("--project" "$2")
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
        -k|--k)
            check_usage "-k/--k" $2
            pre_opts+=("--k" "$2")
            shift 2;;
        --minkf)
            check_usage "--minkf" $2
            pre_opts+=("--minkf" "$2")
            shift 2;;
        --maxkf)
            check_usage "--maxkf" $2
            pre_opts+=("--maxkf" "$2")
            shift 2;;
        --correlation-thresh)
            check_usage "--correlation-thresh" $2 
            pre_opts+=("--correlation-thresh" "$2")
            shift 2;;
        --postgres)
            check_usage "--postgres" $2
            psl_opts+=("--postgres" "$2")
            shift 2;;
        -p|--param)
            pre_opts+=("--param")
            shift;;
        --truth)
            check_usage "--truth" $2
            pre_opts+=("--truth" "$2")
            psl_opts+=("--weight_learning")
            shift 2;;
        --baseline)
            check_usage "--baseline" $2
            pre_opts+=("--baseline" "$2")
            shift 2;;
        --) shift;;
        *) echo "Invalid argument: $1" ; exit 1;;
    esac
done
pre_opts+=("--threads" "$(get_cpus)")
pre_opts+=("--mem" "$(get_mem)")
psl_opts+=("--mem" "$(get_mem)")
pre_opts=("${proj[@]}" "${pre_opts[@]}")
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
if [[ -r "$OPATH/scored_kmers.txt" ]] \
&& [[ -r "$OPATH/scored_kmers.fsa" ]] \
&& [[ -r "$OPATH/KMERPHENO.txt" ]]; then
    postprocessed=1
    echo -e "$Y postprocessed"
else
    echo -e "$N postprocessed"
fi
# check for preprocessed files
if [[ -r "$PPATH/contains_obs.txt" ]] \
&& [[ -r "$PPATH/block_obs.txt" ]] \
&& [[ -r "$PPATH/kmerPheno_target.txt" ]] \
&& [[ -r "$PPATH/samplePheno_obs.txt" ]] \
&& [[ -r "$PPATH/similarPheno_obs.txt" ]] \
&& [[ -r "$PPATH/similarSample_obs.txt" ]] \
&& [[ -r "$PPATH/dissimilarSample_obs.txt" ]] \
&& [[ -r "$PPATH/pheno_int_map.pkl" ]] \
&& [[ -r "$PPATH/sample_int_map.pkl" ]] \
&& [[ -r "$PPATH/kmer_int_map.pkl" ]] \
&& [[ -r "$PPATH/kmer_sample_map.txt" ]] \
&& [[ -r "$PPATH/kmer_pheno_map.txt" ]] \
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
    if ! [[ -r "$OPATH/KMERPHENO.txt" ]]; then
        run_psl && run_postprocess
    else
        run_postprocess
    fi
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
