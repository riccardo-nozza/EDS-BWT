#!/usr/bin/env python3

import argparse
import random
import re


def load_sequences(file_path):
    sequences = []

    with open(file_path) as f:
        content = f.read()

    # trova tutto ciò che sta tra {}
    blocks = re.findall(r'\{([^}]*)\}', content)

    for block in blocks:
        parts = block.split(",")

        for p in parts:
            seq = p.strip()

            if seq == "":
                continue

            if set(seq) <= set("ACGTN"):
                sequences.append(seq)

    return sequences


def extract_patterns(sequences, length):
    patterns = []

    for seq in sequences:
        if len(seq) < length:
            continue

        for i in range(len(seq) - length + 1):
            patterns.append(seq[i:i+length])

    return patterns


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("input_file", help="input file")
    parser.add_argument("output_file", help="output file")
    parser.add_argument("-l", "--length", type=int, required=True,
                        help="pattern length")
    parser.add_argument("-n", "--num", type=int, default=90,
                        help="number of patterns (default 90)")

    args = parser.parse_args()

    sequences = load_sequences(args.input_file)

    patterns = extract_patterns(sequences, args.length)

    if len(patterns) < args.num:
        raise ValueError("Not enough patterns available")

    selected = random.sample(patterns, args.num)

    with open(args.output_file, "w") as f:
        for i, p in enumerate(selected):
            f.write(p)
            if i != len(selected) - 1:
                f.write("\n")


if __name__ == "__main__":
    main()