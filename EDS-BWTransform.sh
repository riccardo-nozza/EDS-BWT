#!/bin/bash

#eds_to_ebwt.sh NAMEFILE PATHTOOL

GSUFPATH="../../gsufsort"
BUILDPATH="./build"
NAMEFILE=$1
OUTPUT=$2

$BUILDPATH/eds_to_fasta $NAMEFILE.eds $OUTPUT

$GSUFPATH/gsufsort $OUTPUT.fasta --da --bwt --output $OUTPUT.fasta

$BUILDPATH/da_to_everything $OUTPUT

echo "File "$NAMEFILE" done."
