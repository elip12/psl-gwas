#!/bin/bash

INFILE=$1
OUTFILE='data/intermediate/samples.fa'

touch $OUTFILE
cat $INFILE | while read line 
do
    FNAME=`echo "$line" | cut -f 2`
    if [[ -h $FNAME ]]; then
        cat `readlink $FNAME` >> $OUTFILE
    elif [[ -r $FNAME ]]; then
        cat $FNAME >> $OUTFILE
    fi
done

