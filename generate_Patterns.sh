#!/bin/bash

DATASET_DIR=$1
OUTPUT_DIR=$2
SCRIPT=generate_Patterns_from_file.py

for file in "$DATASET_DIR"/*.eds; do

    base=$(basename "$file" .eds)

    for length in 8 16 32; do

        output="${OUTPUT_DIR}/patterns_${base}_${length}.txt"

        python3 "$SCRIPT" "$file" "$output" -l "$length" -n 90

        echo "Generated $output"

    done

done