#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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


def main(ipr_types_file, infile, outfile):
    ipr_types_and_descriptions = load_ipr(ipr_types_file)
    with open(infile, "r") as file_in, open(outfile, "w") as file_out:
        fasta_flag = False
        for line in file_in:
            if not fasta_flag:
                if line.startswith("##FASTA"):
                    fasta_flag = True
                    file_out.write(line)
                elif not line.startswith("#"):
                    contig, tool, feature, start, end, blank1, strand, blank2, col9 = (
                        line.strip().split("\t")
                    )
                    if feature == "CDS":
                        attributes_dict = dict(
                            re.split(r"(?<!\\)=", item)
                            for item in re.split(r"(?<!\\);", col9)
                        )
                        matching_key = next(
                            (
                                key
                                for key in attributes_dict
                                if key.lower() == "interpro"
                            ),
                            None,
                        )
                        if matching_key:
                            # add descriptions to the InterPro terms
                            attributes_dict[matching_key] = add_ipr_descriptions(
                                attributes_dict[matching_key],
                                ipr_types_and_descriptions,
                            )
                        col9_updated = update_col9(attributes_dict)
                        file_out.write(
                            "\t".join(
                                [
                                    contig,
                                    tool,
                                    feature,
                                    start,
                                    end,
                                    blank1,
                                    strand,
                                    blank2,
                                    col9_updated,
                                ]
                            )
                            + "\n"
                        )
                    else:
                        file_out.write(line)
                else:
                    file_out.write(line)
            else:
                file_out.write(line)


def update_col9(attributes_dict):
    return ";".join([f"{key}={value}" for key, value in attributes_dict.items()])


def add_ipr_descriptions(ipr_string, ipr_types_and_descriptions):
    new_ipr_list = list()
    ipr_terms = ipr_string.split(",")
    for ipr_term in ipr_terms:
        if ipr_term in ipr_types_and_descriptions:
            ipr_extended = "{}: {} [{}]".format(
                ipr_term,
                ipr_types_and_descriptions[ipr_term]["desc"],
                ipr_types_and_descriptions[ipr_term]["type"],
            )
            new_ipr_list.append(ipr_extended)
        else:
            new_ipr_list.append(ipr_term)
    return ",".join(new_ipr_list)


def load_ipr(ipr_types_file):
    ipr_types = dict()
    short_types = {
        "Active_site": "S",
        "Binding_site": "S",
        "Conserved_site": "S",
        "Domain": "D",
        "Family": "F",
        "Homologous_superfamily": "H",
        "PTM": "S",
        "Repeat": "R",
    }
    with open(ipr_types_file, "r") as file_in:
        for line in file_in:
            if line.startswith("IPR"):
                acc, ipr_type, desc = line.strip().split("\t")
                ipr_types.setdefault(acc, {}).update(
                    {"type": short_types[ipr_type], "desc": desc}
                )
    return ipr_types


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes an annotated GFF file and adds in descriptions of InterPro terms."
        )
    )
    parser.add_argument(
        "--ipr-entries",
        required=True,
        help="The path to the entries.list file from InterPro.",
    )
    parser.add_argument(
        "-i",
        dest="infile",
        required=True,
        help="The path to the input GFF with all annotations in place, including UniFIRE.",
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
        args.ipr_entries,
        args.infile,
        args.outfile,
    )
