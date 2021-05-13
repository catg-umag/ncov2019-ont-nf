#!/usr/bin/env python3
import argparse
import re

import pandas as pd


def main():
    args = parse_arguments()

    clades = pd.read_csv(args.clades)

    clade_assignments = {}
    for vcf in args.inputs:
        df = load_vcf(vcf)
        name = re.sub(r"^.*/|\..*$", "", vcf)
        clade_assignments[name] = get_clade(df, clades)

    # write output
    with open(args.output, "w") as f:
        f.write("sample,clade\n")
        f.write("\n".join(','.join(x) for x in clade_assignments.items()))


def get_clade(sample_variants: pd.DataFrame, clades: pd.DataFrame) -> str:
    # special case: wildtype (should have . as ref in all its "mutations")
    # so, if we don't find any mutation in the same positions, we assume it's wildtype
    wt = clades.query("ref == '.'")
    wt_name = wt.clade.unique()[0]
    if len(wt.merge(sample_variants, on=["pos"])) == 0:
        return wt_name

    # count how many mutations each clade has
    clade_nmutations = (
        clades.query("clade != @wt_name")
        .groupby("clade")
        .size()
        .to_frame("n")
        .reset_index()
    )

    selected_clade = ("None", 0)
    for t in clade_nmutations.itertuples():
        name = t.clade
        matched = len(sample_variants.merge(clades.query("clade == @name")))
        # if we find all the mutations of this clade and the number of mutations
        # is higher than the current selected clade, choose this clade
        if matched == t.n and t.n > selected_clade[1]:
            selected_clade = (name, t.n)

    return selected_clade[0]


def load_vcf(path: str) -> pd.DataFrame:
    df = pd.read_csv(
        path, comment="#", sep="\t", usecols=[1, 3, 4], names=["pos", "ref", "alt"]
    )

    return df


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Assigns clades from VCF using a list with each clade and its mutations"
    )
    parser.add_argument(
        "--clade-mutation-list",
        "-m",
        required=True,
        dest="clades",
        help="CSV file containing list with each clade and its defining mutations,"
        + " cols: (clade,pos,ref,alt)",
    )
    parser.add_argument(
        "--output",
        "-o",
        required=True,
        help="ouput CSV containing each sample with its clade (basename of files will be used as sample)",
    )
    parser.add_argument(
        "--input-vcfs",
        "-i",
        required=True,
        nargs="+",
        dest="inputs",
        help="input VCFs to determine clade",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
