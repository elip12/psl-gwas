#!/bin/bash
###############################################################################
##  postprocess.sh
##  Invokes postprocess.py script, passing all command-line arguments through.
##  run.sh will parse the arguments and only give this script the ones it needs.
###############################################################################
python3 psl-gwas/postprocess.py "$@"

