#!/usr/bin/env python3

import argparse
import csv
import logging

logging.basicConfig(level=logging.INFO)


def main(amr_file, outfile, version):
    with open(amr_file, "r") as file_in, open(outfile, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t")
        writer.writerow(["##gff-version 3"])
        for line in file_in:
            if line.startswith("Protein identifier"):
                continue
            (
                protein_id,
                contig,
                start,
                end,
                strand,
                gene_name,
                seq_name,
                scope,
                element_type,
                element_subtype,
                drug_class,
                drug_subclass,
                _,
            ) = line.strip().split("\t", 12)
            writer.writerow([contig, f"AMRFinderPlus:{version}", "gene", start, end,
                             ".", strand, ".", f"ID={protein_id};gene_name={gene_name};sequence_name={seq_name};"
                                               f"scope={scope};element_type={element_type};"
                                               f"element_subtype={element_subtype};drug_class={drug_class};"
                                               f"drug_subclass={drug_subclass}"])


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes AMRFinderPlus output and parses it to create a standalone GFF."
        )
    )
    parser.add_argument(
        "-i",
        dest="amr_in",
        required=True,
        help="Path to the AMRFinderPlus result file.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help=(
            "Path to the output file."
        ),
    )
    parser.add_argument(
        "-v",
        dest="version",
        required=True,
        help=(
            "AMRFinderPlus version."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.amr_in,
        args.outfile,
        args.version
    )
