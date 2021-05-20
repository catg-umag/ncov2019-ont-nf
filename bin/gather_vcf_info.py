#!/usr/bin/env python3
import argparse
import re
from typing import List


def main():
    args = parse_arguments()

    variants = []
    for filename in args.inputs:
        variants += process_vcf_file(filename)

    with open(args.output, "w") as f:
        f.write("sample,pos,ref,alt,variant_type,gene,nt_change,aa_change\n")
        for v in variants:
            f.write(",".join(v) + "\n")


def process_vcf_file(filename: str) -> List[List[str]]:
    """
    Read VCF file and extract a number of relevant files for each variant
    """
    ann_fields = []
    variants = []
    name = re.sub(r"^.*/|\..*$", "", filename)
    with open(filename) as f:
        for row in f.readlines():
            if row.startswith("##INFO=<ID=ANN"):
                ann_fields = (
                    re.search(r"(?<=')[A-Za-z|/_. ]+(?=')", row).group().split(" | ")
                )
            elif not row.startswith("#"):
                variants.append([name] + process_vcf_row(row, ann_fields))

    return variants
    

def process_vcf_row(row: str, ann_fields: List[str]) -> List[str]:
    """
    Process a VCF row extracting relevant columns

    :param row: raw VCF row
    :param ann_fields: list of ANN fields
    :return list with selected fields
    """
    ann_selected_fields = ("Annotation", "Gene_Name", "HGVS.c", "HGVS.p")
    vcf_fields = row.split("\t")
    ann_values = re.search(r"(?<=ANN=)([^;]+)", vcf_fields[7]).group().split("|")

    values = [x for i, x in enumerate(vcf_fields, start=1) if i in (2, 4, 5)]
    values += [v for k, v in zip(ann_fields, ann_values) if k in ann_selected_fields]

    return values


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Collects info from VCF and outputs it in a CSV file"
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
        help="input VCFs",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
