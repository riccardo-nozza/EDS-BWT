#!/bin/bash

#eds_to_ebwt.sh NAMEFILE PATHTOOL

NAMEFILE=$1
PATHTOOL=$2

./eds_to_fasta $NAMEFILE.eds $NAMEFILE

"$PATHTOOL"/gsufsort $NAMEFILE.fasta --da --bwt --output $NAMEFILE.fasta

./da_to_everything $NAMEFILE

echo "File "$NAMEFILE" done."
