#!/usr/bin/env python3
import argparse

import pandas as pd


def main():
    args = parse_arguments()

    df_base = pd.read_csv(args.base)
    df_base.drop(
        columns=[x for x in df_base.columns if x.startswith("gisaid_")], inplace=True
    )
    df_lineage = pd.read_csv(args.lineages)
    df_coverage = pd.read_csv(
        args.coverage, sep="\t", usecols=["#rname", "meandepth"]
    ).rename(columns={"#rname": "sample"})
    df_refcoverage = pd.read_csv(args.ref_coverage)

    df = (
        df_base.merge(df_lineage, how="outer")
        .merge(df_coverage, how="outer")
        .merge(df_refcoverage, how="outer")
    ).sort_values("sample")

    df.to_csv(args.output, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merges multiple data (base, lineages, coverage, ref_covered)"
    )

    parser.add_argument(
        "--base-data",
        "-b",
        required=True,
        dest="base",
        help="CSV file with at leas sample and barcode",
    )
    parser.add_argument("--lineages", "-l", required=True, help="Lineages CSV file")
    parser.add_argument(
        "--coverage-stats",
        "-c",
        required=True,
        dest="coverage",
        help="Coverage stats TSV file",
    )
    parser.add_argument(
        "--ref-coverage-stats",
        "-r",
        required=True,
        dest="ref_coverage",
        help="Reference coverage stats CSV file",
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")

    return parser.parse_args()


if __name__ == "__main__":
    main()
