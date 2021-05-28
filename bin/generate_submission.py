#!/usr/bin/env python3
import argparse
import csv
import re
from typing import Dict

from Bio import SeqIO
from openpyxl import load_workbook
from openpyxl.worksheet.worksheet import Worksheet


def main():
    args = parse_arguments()

    sample_data = prepare_sample_info(
        args.base_data, args.sample_summary, args.output_sequences, args.min_coverage
    )

    # load sequences
    process_sequences(args.sequences, args.output_sequences, sample_data)

    # load template
    wb = load_workbook(filename=args.template)
    fill_sheet(wb["Submissions"], sample_data)
    wb.save(filename=args.output)


def fill_sheet(sheet: Worksheet, sample_data: Dict[str, Dict[str, str]]):
    # delete example row
    sheet.delete_rows(3)

    # get column letter for each field
    columns = {col[0].value: col[0].column_letter for col in sheet.columns}

    if "covv_virus_name" in sample_data[list(sample_data.keys())[0]]:
        data = sorted(sample_data.values(), key=lambda x: x["covv_virus_name"])
    else:
        data = [v for k, v in sorted(sample_data.items(), key=lambda x: x[0])]

    for i, row in enumerate(data, start=3):
        for k, v in row.items():
            sheet[f"{columns[k]}{i}"].value = v


def prepare_sample_info(
    base_data: str,
    sample_summary: str,
    output_fasta_filename: str,
    min_coverage_perc: int,
) -> Dict[str, Dict[str, str]]:
    sample_data = {}

    # load coverages
    with open(sample_summary) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if float(row["perc_covered"]) >= min_coverage_perc:
                coverage = round(float(row["meandepth"]))
                sample_data[row["sample"]] = {"covv_coverage": f"{coverage}x"}

    # load base data
    with open(base_data) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["sample"] in sample_data:
                for k, v in row.items():
                    if k.startswith("gisaid_"):
                        sample_data[row["sample"]][k.replace("gisaid_", "")] = v

    # add generic info
    for info in sample_data.values():
        info["fn"] = output_fasta_filename
        info["covv_type"] = "betacoronavirus"
        info["covv_passage"] = "Original"
        info["covv_host"] = "Human"
        info["covv_seq_technology"] = "Nanopore MinION"
        info["covv_assembly_method"] = "ARTIC-nCoV-bioinformaticsSOP-v1.1.0"

    return sample_data


def process_sequences(
    input_filename: str, output_filename: str, sample_data: Dict[str, Dict[str, str]]
):
    """
    Load sequences, filter presence in sample_data, rename and write them

    :param input_filename: current filename with sequences
    :param output_filename: output filename
    :param sample_data: sample data (covv_virus_name will be used if defined)
    """
    # load sequences
    sequences = []
    for s in SeqIO.parse(input_filename, "fasta"):
        name = re.sub("/ARTIC.*", "", s.id)

        if name not in sample_data:
            continue

        if "covv_virus_name" in sample_data[name]:
            s.id = s.description = sample_data[name]["covv_virus_name"]
        else:
            s.id = s.description = name

        sequences.append(s)

    sequences = sorted(sequences, key=lambda x: x.id)

    SeqIO.write(sequences, output_filename, "fasta")


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Prepares GISAID submission filling sheet "
        + "and generating FASTA file with sequences"
    )

    parser.add_argument(
        "--submission-template",
        "-t",
        required=True,
        dest="template",
        help="EpiCoV submission template (.xlsx)",
    )
    parser.add_argument(
        "--sample-summary",
        "-s",
        required=True,
        dest="sample_summary",
        help="sample summary (.csv)",
    )
    parser.add_argument(
        "--base-data",
        "-b",
        required=True,
        dest="base_data",
        help="base sample data (with GISAID fields) (.csv)",
    )
    parser.add_argument(
        "--input-sequences",
        "-i",
        required=True,
        dest="sequences",
        help="consensus for each sample (.fasta)",
    )
    parser.add_argument(
        "--required-ref-coverage",
        "-r",
        default=95,
        type=int,
        dest="min_coverage",
        help="required reference coverage (as percentage) to call a valid assembly",
    )
    parser.add_argument(
        "--output-excel",
        "-o",
        required=True,
        dest="output",
        help="Output filled sheet (.xlsx)",
    )
    parser.add_argument(
        "--output-sequences",
        "-O",
        required=True,
        dest="output_sequences",
        help="Output sequences (.fasta)",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
