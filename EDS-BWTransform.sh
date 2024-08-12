#!/bin/bash

#eds_to_ebwt.sh NAMEFILE PATHTOOL

GSUFPATH="gsufsort"
NAMEFILE=$1
OUTPUT=$2

./eds_to_fasta $NAMEFILE.eds $OUTPUT

$GSUFPATH/gsufsort $NAMEFILE.fasta --da --bwt --output $OUTPUT.fasta

./da_to_everything $OUTPUT

echo "File "$NAMEFILE" done."
