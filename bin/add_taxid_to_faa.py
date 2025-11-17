#!/usr/bin/env python3

# Copyright 2023-2025 EMBL - European Bioinformatics Institute
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
import logging
import sys

logging.basicConfig(level=logging.INFO)


def main(infile, taxid, outfile):
    if not taxid.isdigit():
        sys.exit(
            f"Taxid must consist of digits only. Taxid {taxid} is not valid. Exiting."
        )

    with open(outfile, "w") as file_out, open(infile) as file_in:
        for line in file_in:
            if line.startswith(">"):
                formatted_line = reformat_line(line, taxid)
                file_out.write(formatted_line)
            else:
                file_out.write(line)


def reformat_line(line, taxid):
    line = line.lstrip(">").strip()
    protein_id, description = line.split(maxsplit=1)
    description = (
        description.replace('"', "").replace("'", "").replace("‘", "").replace("’", "")
    )
    formatted_line = f">tr|{protein_id}|{description} OX={taxid}\n"
    return formatted_line


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script reformats the fasta faa file to prepare it for UniFIRE."
        )
    )
    parser.add_argument(
        "-i",
        dest="infile",
        required=True,
        help="Input protein fasta file.",
    )
    parser.add_argument(
        "-t",
        dest="taxid",
        required=True,
        help="NCBI taxid for the species or the lowest taxonomy rank known for the genome.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="Path to the file where the output will be saved to.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.infile,
        args.taxid,
        args.outfile,
    )
