#!/usr/bin/env python3
import argparse
import re

from Bio import SeqIO


GENOME_SIZE = 29903


def main():
    args = parse_arguments()

    stats = []
    for seq in SeqIO.parse(args.input, "fasta"):
        name = re.sub(r"/ARTIC.*$", "", seq.id)
        n_counts = len(re.findall("N", str(seq.seq)))
        stats.append([name, n_counts])

    stats = sorted(stats, key=lambda x: x[0])

    with open(args.output, "w") as f:
        f.write("sample,bases_covered,perc_covered\n")
        for s in stats:
            covered = round(100 - (s[1] / GENOME_SIZE * 100), 2)
            f.write(f"{s[0]},{GENOME_SIZE - s[1]},{covered}\n")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Collect statistics regarding the completion of the assembly "
        + "compared to the reference (reference coverage)"
    )

    parser.add_argument(
        "--input", "-i", required=True, help="FASTA file containing all consensuses"
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")

    return parser.parse_args()


if __name__ == "__main__":
    main()
