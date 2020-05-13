#!/bin/bash
###############################################################################
##  run_psl.sh
##  Invokes PSL. 
##  run.sh will parse the arguments and only give this script the ones it needs.
###############################################################################
readonly PSL_VERSION='2.3.0-SNAPSHOT'
readonly JAR_PATH="./psl-cli-${PSL_VERSION}.jar"
readonly ADDITIONAL_PSL_OPTIONS='-D log4j.threshold=TRACE --int-ids -D eval.closetruth=false -D runtimestats.collect=true --skipAtomCommit'
readonly ADDITIONAL_LEARN_OPTIONS='--learn ContinuousRandomGridSearch -D weightlearning.evaluator=ContinuousEvaluator -D continuousrandomgridsearch.maxlocations=250'
readonly ADDITIONAL_EVAL_OPTIONS='--infer SGDStreamingInference --eval org.linqs.psl.evaluation.statistics.ContinuousEvaluator' 

BASE_NAME=''
MEM=''
runweightlearning=0

function main() {
    trap exit SIGINT

    ARGS=()
    while (( "$#" )) ; do
        case "$1" in
            --project) BASE_NAME="$2" ; shift 2;;
            --mem) MEM="$2G" ; shift 2;;
            --weight_learning) runweightlearning=1 ; shift;;
            *) ARGS+=("$1") ; shift ;;
        esac
    done
    set -- ${ARGS[@]}
    if [[ -z BASE_NAME ]]; then
        echo "No project name given. Aborting."
        exit 1
    fi

    # Make sure we can run PSL.
    check_requirements
    fetch_psl

    # Run PSL
    if [[ runweightlearning -eq 1 ]]; then
        runWeightLearning "$@"
    fi
    runEvaluation "$@"
}


function runWeightLearning() {
    echo "Running PSL Weight Learning"

    java -Xmx350G -Xms350G -jar "${JAR_PATH}" --model "${BASE_NAME}/gwas.psl" --data "${BASE_NAME}/gwas.data" ${ADDITIONAL_LEARN_OPTIONS} ${ADDITIONAL_PSL_OPTIONS} "$@"
    if [[ "$?" -ne 0 ]]; then
        echo 'ERROR: Failed to run weight learning'
        exit 60
    fi
}

function runEvaluation() {
    echo "Running PSL Inference"

    java -Xmx350G -Xms350G -jar "${JAR_PATH}" --model "${BASE_NAME}/gwas.psl" --data "${BASE_NAME}/gwas.data" --output ${BASE_NAME}/data/postprocessed ${ADDITIONAL_EVAL_OPTIONS} ${ADDITIONAL_PSL_OPTIONS} "$@"
    if [[ "$?" -ne 0 ]]; then
        echo 'ERROR: Failed to run infernce'
        exit 70
    fi
}

function check_requirements() {
    local hasWget
    local hasCurl

    type wget > /dev/null 2> /dev/null
    hasWget=$?

    type curl > /dev/null 2> /dev/null
    hasCurl=$?

    if [[ "${hasWget}" -ne 0 ]] && [[ "${hasCurl}" -ne 0 ]]; then
        echo 'ERROR: wget or curl required to download dataset'
        exit 10
    fi

    type java > /dev/null 2> /dev/null
    if [[ "$?" -ne 0 ]]; then
        echo 'ERROR: java required to run project'
        exit 13
    fi
}

function get_fetch_command() {
    type curl > /dev/null 2> /dev/null
    if [[ "$?" -eq 0 ]]; then
        echo "curl -o"
        return
    fi

    type wget > /dev/null 2> /dev/null
    if [[ "$?" -eq 0 ]]; then
        echo "wget -O"
        return
    fi

    echo 'ERROR: wget or curl not found'
    exit 20
}

function fetch_file() {
    local url=$1
    local path=$2
    local name=$3

    if [[ -e "${path}" ]]; then
        echo "${name} file found cached, skipping download."
        return
    fi

    echo "Downloading ${name} file located at: '${url}'."
    `get_fetch_command` "${path}" "${url}"
    if [[ "$?" -ne 0 ]]; then
        echo "ERROR: Failed to download ${name} file"
        exit 30
    fi
}

# Fetch the jar from a remote or local location and put it in this directory.
# Snapshots are fetched from the local maven repo and other builds are fetched remotely.
function fetch_psl() {
    if [[ $PSL_VERSION == *'SNAPSHOT'* ]]; then
        local snapshotJARPath="$HOME/.m2/repository/org/linqs/psl-cli/${PSL_VERSION}/psl-cli-${PSL_VERSION}.jar"
        cp "${snapshotJARPath}" "${JAR_PATH}"
    else
        local remoteJARURL="https://repo1.maven.org/maven2/org/linqs/psl-cli/${PSL_VERSION}/psl-cli-${PSL_VERSION}.jar"
        fetch_file "${remoteJARURL}" "${JAR_PATH}" 'psl-jar'
    fi
}

main "$@"
