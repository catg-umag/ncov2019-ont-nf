#!/usr/bin/env python3
import argparse
import re

from Bio import SeqIO


def main():
    args = parse_arguments()

    stats = []
    max_length = 0
    for seq in SeqIO.parse(args.input, "fasta"):
        name = re.sub(r"/ARTIC.*$", "", seq.id)
        n_counts = len(re.findall("N", str(seq.seq)))
        if len(seq.seq) > max_length:
            max_length = len(seq.seq)
        stats.append([name, n_counts])
    
    stats = sorted(stats, key = lambda x: int(x[0]))

    with open(args.output, "w") as f:
        f.write("sample,n_quantity,perc_covered\n")
        for s in stats:
            covered = round(100 - (s[1] / max_length * 100), 2)
            f.write(f"{s[0]},{s[1]},{covered}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Counts how many 'N' are found in each sequecing, and generate stats"
    )

    parser.add_argument(
        "--input", "-i", required=True, help="FASTA file containing all consensuses"
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")

    return parser.parse_args()


if __name__ == "__main__":
    main()
