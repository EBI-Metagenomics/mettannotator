#!/usr/bin/env python3

# Copyright 2023-2024 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import argparse
import csv
import logging

logging.basicConfig(level=logging.INFO)


def main(amr_file, outfile, version):
    with open(amr_file) as file_in, open(outfile, "w") as file_out:
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
            writer.writerow(
                [
                    contig,
                    f"AMRFinderPlus:{version}",
                    "gene",
                    start,
                    end,
                    ".",
                    strand,
                    ".",
                    f"ID={protein_id};gene_name={gene_name};sequence_name={seq_name};"
                    f"scope={scope};element_type={element_type};"
                    f"element_subtype={element_subtype};drug_class={drug_class};"
                    f"drug_subclass={drug_subclass}",
                ]
            )


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
        help=("Path to the output file."),
    )
    parser.add_argument(
        "-v",
        dest="version",
        required=True,
        help=("AMRFinderPlus version."),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args.amr_in, args.outfile, args.version)
