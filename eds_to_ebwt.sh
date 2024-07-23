#!/bin/bash

#eds_to_ebwt.sh NAMEFILE PATHTOOL

NAMEFILE=$1
PATHTOOL=$2

./eds_to_fasta $NAMEFILE.eds $NAMEFILE

gsufsort $NAMEFILE.fasta --da --bwt --output $NAMEFILE.fasta

chmod 777 $NAMEFILE.empty.sh
$NAMEFILE.empty.sh
rm $NAMEFILE.empty.sh

./da_to_everything $NAMEFILE

chmod 777 $NAMEFILE.ebwt
chmod 777 $NAMEFILE.bitvector
chmod 777 $NAMEFILE*.aux
rm $NAMEFILE.fasta
rm $NAMEFILE.fasta.4.da
rm $NAMEFILE.fasta.bwt

echo "File "$NAMEFILE" done."
