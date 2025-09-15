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


def process_file(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            cols = line.rstrip("\n").split("\t")
            if cols:  # make sure line isn’t empty
                # value in cols[0] will look like tr|GCA_031506065_00606|hypothetical
                # we want to only keep the actual accession (GCA_031506065_00606)
                parts = cols[0].split("|")
                if len(parts) >= 2:
                    cols[0] = parts[1]  # take the part between first and second "|"
                    # if there are no pipes in the string, the value is unchanged
            outfile.write("\t".join(cols) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Clean protein accessions in IPS output files.")
    parser.add_argument("-i", "--input", required=True, help="Input tab-delimited IPS file")
    parser.add_argument("-o", "--output", required=True, help="Output file name")
    args = parser.parse_args()

    process_file(args.input, args.output)


if __name__ == "__main__":
    main()
