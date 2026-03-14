#!/usr/bin/env python3

import argparse
import random


def load_sequences(msa_file):

    sequences = []
    seq = []

    with open(msa_file) as f:
        for line in f:

            if line.startswith(">"):
                if seq:
                    sequences.append("".join(seq))
                    seq = []
                continue

            cleaned = line.strip().replace("-", "")
            seq.append(cleaned)

    if seq:
        sequences.append("".join(seq))

    return sequences


def extract_patterns(sequences, length):

    patterns = []

    for s in sequences:

        if len(s) < length:
            continue

        for i in range(len(s) - length + 1):
            patterns.append(s[i:i+length])

    return patterns


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("input")
    parser.add_argument("output")
    parser.add_argument("-l", "--length", type=int, required=True)
    parser.add_argument("-n", "--num", type=int, default=90)

    args = parser.parse_args()

    sequences = load_sequences(args.input)

    patterns = extract_patterns(sequences, args.length)

    selected = random.sample(patterns, args.num)

    with open(args.output, "w") as f:
        for i, p in enumerate(selected):

            f.write(p)

            if i != len(selected)-1:
                f.write("\n")


if __name__ == "__main__":
    main()