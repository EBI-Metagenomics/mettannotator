#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2023 EMBL - European Bioinformatics Institute
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


def split_cds_to_gene_exon_mrna(entry):
    """Function to split CDS entries into gene, exon, and mRNA entries"""
    if len(entry) > 8 and entry[2] == "CDS":
        # replicate the line
        gene_entry = entry[:]
        exon_entry = entry[:]
        mrna_entry = entry[:]
        cds_entry = entry[:]

        # Update the feature type
        gene_entry[2] = "gene"
        exon_entry[2] = "exon"
        mrna_entry[2] = "mRNA"
        cds_entry[2] = "CDS"

        # tweak IDs and add Parent attributes
        attributes = entry[8]
        attributes_dict = dict(item.split("=") for item in attributes.split(";"))
        if "ID" in attributes_dict:
            gene_id = attributes_dict["ID"]  # use the cds ID as the gene ID

        gene_entry[8] = edit_id_add_parent(attributes_dict, gene_id)
        mrna_entry[8] = edit_id_add_parent(
            attributes_dict, f"transcript:{gene_id}", gene_id
        )
        exon_entry[8] = edit_id_add_parent(
            attributes_dict, f"exon:{gene_id}", f"transcript:{gene_id}"
        )
        cds_entry[8] = edit_id_add_parent(
            attributes_dict, f"CDS:{gene_id}", f"transcript:{gene_id}"
        )

        return [gene_entry, mrna_entry, exon_entry, cds_entry]
    else:
        return [entry]  # comments, ncrna lines, sequences, possibly also invalid gff!


def edit_id_add_parent(attributes_dict, id, parent=None):
    """Function to edit the ID and add Parent attribute"""
    attributes_dict["ID"] = id
    if parent is not None:
        attributes_dict["Parent"] = parent
    return ";".join(f"{key}={value}" for key, value in attributes_dict.items())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Split CDS entries in a GFF file into gene, exon, and mRNA entries"
    )
    parser.add_argument("-i", dest="input_file", required=True, help="Input GFF file")
    parser.add_argument("-o", dest="output_file", required=True, help="Output GFF file")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    # Read the input GFF file and process the entries
    output_entries = []
    with open(input_file, "r") as infile:
        for line in infile:
            entry = line.strip().split("\t")
            new_entries = split_cds_to_gene_exon_mrna(entry)
            output_entries.extend(new_entries)

    # Write the modified GFF entries to the output file
    with open(output_file, "w") as outfile:
        for entry in output_entries:
            outfile.write("\t".join(entry) + "\n")
