#!/usr/bin/env python3

# Copyright 2024 EMBL - European Bioinformatics Institute
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
import re


def main(input_file, output_file):
    with open(input_file, "r") as file_in, open(output_file, "w") as file_out:
        for line in file_in:
            if line.startswith("#") or not line.strip():
                # Write header or empty lines as is
                file_out.write(line)
                continue

            columns = line.strip().split("\t")
            if len(columns) < 9:
                file_out.write(line)
                continue

            if columns[2] != "tRNA":
                file_out.write(line)
                continue

            # Parse the attributes (column 9)
            attributes = columns[8].rstrip(";")
            attributes_dict = dict(
                re.split(r"(?<!\\)=", item)
                for item in re.split(r"(?<!\\);", attributes)
            )

            # Add locus_tag based on the ID field
            if "ID" in attributes_dict:
                locus_tag = attributes_dict["ID"]
                attributes = attributes + f";locus_tag={locus_tag}"

            columns[8] = attributes
            file_out.write("\t".join(columns) + "\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description=("The script adds locus tags to the 9th column.")
    )
    parser.add_argument(
        "-i",
        dest="infile",
        required=True,
        help="The path to the input GFF.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="Path to the output file where the result will be saved.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.infile,
        args.outfile,
    )
