#!/usr/bin/env python3
import argparse

import pandas as pd


def main():
    args = parse_arguments()

    df_pangolin = (
        pd.read_csv(
            args.pangolin,
            usecols=["taxon", "lineage", "scorpio_call"],
            dtype={"taxon": str, "lineage": str, "scorpio_call": str},
        )
        .rename(
            columns={
                "taxon": "sample",
                "lineage": "lineage_pango",
                "scorpio_call": "variant",
            }
        )
        .assign(
            sample=lambda x: x["sample"].str.replace(r"/ARTIC.*$", "", regex=True),
            variant=lambda x: x["variant"].replace(r" \(.*", "", regex=True),
        )
    )
    df_nextstrain = (
        pd.read_csv(
            args.nextstrain,
            usecols=["seqName", "clade"],
            sep="\t",
            dtype={"seqName": str, "clade": str},
        )
        .rename(columns={"seqName": "sample", "clade": "clade_nextstrain"})
        .assign(
            sample=lambda x: x["sample"].str.replace(r"/ARTIC.*$", "", regex=True),
            clade_nextstrain=lambda x: x["clade_nextstrain"].str.replace(
                " .*$", "", regex=True
            ),
        )
    )

    df = df_pangolin.merge(df_nextstrain).sort_values(["sample"])

    df.to_csv(args.output, index=False)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Merges lineage info from GISAID, Pangolin and Nexstrain in a single file"
    )

    parser.add_argument(
        "--pangolin", "-p", required=True, help="Pangolin lineages CSV file"
    )
    parser.add_argument(
        "--nextstrain", "-n", required=True, help="Nextrain lineages CSV file"
    )
    parser.add_argument("--output", "-o", required=True, help="output CSV file")

    return parser.parse_args()


if __name__ == "__main__":
    main()
