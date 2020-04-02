#!/bin/bash
###############################################################################
# startproject.sh
# This script initializes a PSL-GWAS project. The user specifies a name
# for the project, which will have separate data and logs from any other
# project. The preprocessing, association testing, and postprocessing
# code gets shared between all projects.
###############################################################################

# copies .parameters.yaml to project parameters.yaml for user to edit,
# removing do not edit message at top
project_params() {
    tail -n +3 .parameters.yaml > $1
    chmod o+w $1
}

gwas_data() {
    cp psl-gwas/psl/gwas.psl $1/gwas.psl
    sed "s/data/$1\/data/g" psl-gwas/psl/gwas.data > $1/gwas.data
    chmod o+w $1/gwas.psl $1/gwas.data
}

# ensure user inputs name
if [[ "$#" -ne 1 ]]; then
    echo "Please specify a project name."
    exit 1
fi
echo "Initializing project $1" 
echo
# check for python3 version greater than 3.6.0
echo "Checking python3 version..."
version=$(/usr/bin/env python3 -V 2>&1 | grep -Go "[0-9]\.[0-9]\.[0-9]")
intversion=$(echo "$version" | sed -e "s/\.//g")
if [[ -z "$intversion" ]]; then
    echo "No python3 distribution found"
    exit 1
elif [[ "$intversion" -lt 360 ]]; then
    echo "python3 version must be >= 3.6.0"
    exit 1
fi
echo "Found python3 distribution: $version."
echo
# advise user to create a python virtual environment and prompt to continue
echo "If you have not created and activated a python3 virtual environment,"
echo "now is the time to do so."
read -p "Continue? [y/n] " yn
case $yn in
    [Yy] ) ;;
    [Nn] ) exit;;
    * ) echo "Please answer y or n.";;
esac
echo
# install/update python modules
echo "Updating python3 modules..."
python3 -m pip install -U pip --quiet
python3 -m pip install -r requirements.txt --quiet
echo 'python3 modules up to date.'
echo
# create project data directory. All data files PSL-GWAS generates will
# be stored in these directories
echo "Creating project data directory and subdirectories..."
mkdir $1
mkdir $1/data
mkdir $1/data/raw
mkdir $1/data/preprocessed
mkdir $1/data/postprocessed
echo "Successfully created project data directory and subdirectories."
echo
echo "Creating project logs directory..."
mkdir $1/logs
echo "Successfully created project logs directory."
echo
echo "Creating parameters file..."
project_params "$1/parameters.yaml"
echo "Sucessfully created parameters file."
echo
echo "Creating project psl files..."
gwas_data "$1"
echo "Successfully created project psl files."
echo 'Done.'

