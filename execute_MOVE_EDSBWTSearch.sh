#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <mode:0|1> <input_file_or_directory> <output> <patterns_file> <alpha>" >&2
    exit 1
fi

MODE="$1"
INPUT="$2"
OUTPUT="$3"
PATTERNSFILE="$4"
ALPHA="$5"
BUILDPATH="./build"

run_pipeline() {
    local NAMEFILE="$1"
    local OUTFILE="$2"

    echo "----------------------------------------"
    echo "Processing: $NAMEFILE"
    echo "----------------------------------------"

    echo "[1/3] Running EDS-BWTransform..."
    ./EDS-BWTransform.sh "$NAMEFILE" "$OUTFILE"

    echo "[2/3] Building M_LF structure..."
    "$BUILDPATH/build_MLF" "$OUTFILE" "$ALPHA"

    echo "[3/3] Running MOVE_EDSBWTSearch..."
    "$BUILDPATH/MOVE_EDSBWTSearch" "$OUTFILE" "$PATTERNSFILE"

    echo "Completed: $NAMEFILE"
}

echo "[execute_MOVE_EDSBWTSearch] Starting..."
echo "Mode  : $MODE"
echo "Alpha : $ALPHA"

# =========================
# MODE = Single file
# =========================
if [ "$MODE" -eq 0 ]; then

    if [ ! -f "$INPUT.eds" ]; then
        echo "Error: Input file does not exist." >&2
        exit 1
    fi

    run_pipeline "$INPUT" "$OUTPUT"

# =========================
# MODE 1 = All files in directory
# =========================
elif [ "$MODE" -eq 1 ]; then

    if [ ! -d "$INPUT" ]; then
        echo "Error: Input directory does not exist." >&2
        exit 1
    fi

    for file in "$INPUT"*; do
        filename=$(basename "$file")
        filename_noext="${filename%.eds}"
        outfile="${OUTPUT}_${filename_noext}"
        run_pipeline "$INPUT$filename_noext" "$outfile"
    done

else
    echo "Error: mode must be 0 (single file) or 1 (directory)" >&2
    exit 1
fi

echo "All done."