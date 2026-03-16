#!/bin/bash


BCRPATH="BCR_LCP_GSA"
GSUFPATH="../../gsufsort"

NAMEFILE=$1
OUTPUT=$2
BCR=$3

if [[ "$(tail -c 1 "$NAMEFILE")" != "}" ]]; then
    echo "ERROR: file .eds must end with }"
    exit 1
fi

./eds_to_fasta $NAMEFILE $OUTPUT

if [ $BCR -eq 1 ]
then
    $BCRPATH/"BCR_LCP_GSA" $OUTPUT.fasta $OUTPUT
	rm $OUTPUT.fasta
	rm $OUTPUT.len
	rm $OUTPUT.info
	./EOFpos_to_everything $OUTPUT
else
$GSUFPATH/"gsufsort" $OUTPUT.fasta --da --bwt --output $OUTPUT
	rm $OUTPUT.fasta
	./da_to_everything $OUTPUT
fi 
#$GSUFPATH/gsufsort $OUTPUT.fasta --da --bwt --output $OUTPUT.fasta
echo "File "$NAMEFILE" done."
