#!/bin/bash

DATASET_DIR=$1
OUTPUT_DIR=$2
SCRIPT=extract_patterns_from_msa.py

for file in "$DATASET_DIR"/*.msa; do

    base=$(basename "$file" .msa)

    for length in 8 16 32 64; do

        output="${OUTPUT_DIR}/patterns_${base}_${length}.txt"

        python3 "$SCRIPT" "$file" "$output" -l "$length" -n 100

        echo "Generated $output"

    done

done