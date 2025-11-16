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

import argparse


def process_file(input_file, output_file):
    with open(input_file, encoding="utf-8") as infile:
        lines = infile.readlines()

    # Process the header line - fields look like "Protein id", "Contig id", etc
    # After replacements below, they are changed to "Protein_id", "Contid_id" (requirement for import into the portal)
    header = lines[0].strip().split("\t")
    header = [field.replace(" ", "_") for field in header]
    header_line = "\t".join(header) + "\n"

    # Write to output file
    with open(output_file, "w", encoding="utf-8") as outfile:
        outfile.write(header_line)
        # Write the rest of the lines unchanged
        outfile.writelines(lines[1:])


def main():
    parser = argparse.ArgumentParser(
        description="Replaces spaces in the header of the AMRFinderPlus TSV output file with underscores."
    )
    parser.add_argument("-i", "--input", required=True, help="Input file path")
    parser.add_argument("-o", "--output", required=True, help="Output file path")

    args = parser.parse_args()

    process_file(args.input, args.output)


if __name__ == "__main__":
    main()
