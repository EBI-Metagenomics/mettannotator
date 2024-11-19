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
import os
import re
import subprocess
import sys


def load_uniprot_descriptions(uniprot_tsv):
    uniprot_descriptions = dict()
    entry_index = None
    protein_name_index = None
    with open(uniprot_tsv) as file_in:
        for line in file_in:
            if line.startswith("Entry"):
                fields = line.strip().split("\t")
                for index, field in enumerate(fields):
                    if field == "Entry":
                        entry_index = index
                    elif field == "Protein names":
                        protein_name_index = index
                if protein_name_index is None or entry_index is None:
                    sys.exit("Could not get uniprot indices to load descriptions")
            else:
                fields = line.strip().split("\t")
                uniprot_descriptions[fields[entry_index]] = escape_reserved_characters(
                    fields[protein_name_index]
                )
    assert len(uniprot_descriptions) > 0, "Uniprot description dictionary is empty."
    return uniprot_descriptions


def escape_reserved_characters(string):
    reserved_characters = [";", "=", "&"]
    for ch in reserved_characters:
        if ch in string:
            if ch == ";":
                string = string.replace(ch, "/")
            else:
                string = string.replace(ch, f"\{ch}")
    return string


def run_blast(pipeline_fasta, uniprot_fasta, blast_output):
    # Create blast database
    subprocess.run(
        ["makeblastdb", "-in", pipeline_fasta, "-parse_seqids", "-dbtype", "prot"]
    )

    # Run blastp
    subprocess.run(
        [
            "blastp",
            "-db",
            pipeline_fasta,
            "-query",
            uniprot_fasta,
            "-out",
            blast_output,
            "-outfmt",
            "6",
            "-evalue",
            "1e-10",
        ]
    )


def extract_best_hits(blast_output):
    best_hits = {}
    with open(blast_output) as file:
        for line in file:
            fields = line.strip().split("\t")
            query_id = fields[0]
            result_id = fields[1]
            percent_id = fields[2]
            overlap = fields[3]
            if query_id not in best_hits or float(percent_id) > float(
                best_hits[query_id][2]
            ):
                best_hits[query_id] = (result_id, percent_id, overlap)
    return best_hits


def create_mapping_file(best_hits, mapping_file):
    with open(mapping_file, "w") as file:
        file.write("Query\tResult\tPercent_ID\tOverlap\n")
        for query_id, values in best_hits.items():
            result_id, percent_id, overlap = values
            file.write(f"{query_id}\t{result_id}\t{percent_id}\t{overlap}\n")


def update_gff_with_mapping(gff_file, mapping_file, output_file, uniprot_descriptions):
    mapping_dict = {}
    with open(mapping_file) as file:
        next(file)  # Skip header
        for line in file:
            query, result, _, _ = line.strip().split("\t")
            mapping_dict[result] = query

    with open(gff_file) as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fields = line.strip().split("\t")
                if (
                    len(fields) >= 9
                    and fields[2] == "CDS"
                    and fields[8].startswith("ID=")
                ):
                    attributes_dict = dict(
                        re.split(r"(?<!\\)=", item)
                        for item in re.split(r"(?<!\\);", fields[8])
                    )
                    identifier = attributes_dict["ID"]
                    if identifier.startswith("CDS:"):
                        identifier = identifier[len("CDS:") :]
                    if identifier in mapping_dict:
                        if "Dbxref" in attributes_dict:
                            attributes_dict[
                                "Dbxref"
                            ] += f",UniProt:{mapping_dict[identifier]}"
                        else:
                            attributes_dict["Dbxref"] = (
                                f"UniProt:{mapping_dict[identifier]}"
                            )
                        if mapping_dict[identifier] in uniprot_descriptions:
                            attributes_dict["uniprot_prot_name"] = uniprot_descriptions[
                                mapping_dict[identifier]
                            ]
                        fields[8] = ";".join(
                            [f"{key}={value}" for key, value in attributes_dict.items()]
                        )
                outfile.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Need to download proteomes from UniProt (Buniformis: UP000004110, Pvulgatus: UP000002861) & \
                    strip away any unwanted bits in the headers, and have the new GFF and protein fasta files. \
                    Blastp runs on the two protein fasta files, with evalue -10. The best hits are extracted from the \
                    results to produce a mapping file between the new genes and the old uniprot identifiers. \
                    This mapping is used to insert Dbxref=UniProt: key-value attributes to the ninth column of CDS \
                    lines in the final output. The species identifier is just used to name files.",
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Species identifier to be used in file names (for example, BU_ATCC8492)",
    )
    parser.add_argument(
        "-i",
        dest="gff_file",
        help="Path to the GFF file generated by the pipeline.",
        required=True,
    )
    parser.add_argument(
        "--faa",
        dest="pipeline_fasta",
        help="Protein FASTA file generated by Prokka.",
        required=True,
    )
    parser.add_argument(
        "--uniprot",
        dest="uniprot_fasta",
        help="UniProt protein FASTA file to get the identifiers from.",
        required=True,
    )
    parser.add_argument(
        "-o",
        dest="output_gff_file",
        help="Path to the output file where the results will be saved to. Default: <species>.withuniprotids.output.gff",
        required=False,
    )
    parser.add_argument(
        "--description",
        dest="include_description",
        action="store_true",
        default=False,
        help="Use this flag to add descriptions of UniProt accessions to the output GFF. Must provide UniProt Input.",
        required=False,
    )
    parser.add_argument(
        "--uniprot-tsv",
        help="Provide a path to UniProt TSV if using the --description flag.",
        required=False,
    )

    args = parser.parse_args()

    species = args.species
    if " " in species:
        sys.exit(f"The species parameter {species} is invalid. Spaces are not allowed")
    if args.include_description:
        if not args.uniprot_tsv:
            sys.exit("Must provide a UniProt TSV if the --description flag is used")
        uniprot_tsv = args.uniprot_tsv
        uniprot_descriptions = load_uniprot_descriptions(uniprot_tsv)
    gff_file = args.gff_file
    uniprot_fasta = args.uniprot_fasta
    pipeline_fasta = args.pipeline_fasta
    blast_output = species + ".blastp.evalue.out"
    mapping_file = species + ".mapping.txt"
    if args.output_gff_file:
        output_gff_file = args.output_gff_file
    else:
        output_gff_file = species + ".withuniprotids.output.gff"

    # Run blast if blast results do not already exist
    if not os.path.exists(blast_output):
        run_blast(pipeline_fasta, uniprot_fasta, blast_output)

    # Extract best hits
    best_hits = extract_best_hits(blast_output)

    # Create mapping file
    create_mapping_file(best_hits, mapping_file)

    # Update GFF with mapping
    if args.include_description:
        update_gff_with_mapping(
            gff_file, mapping_file, output_gff_file, uniprot_descriptions
        )
    else:
        update_gff_with_mapping(gff_file, mapping_file, output_gff_file, dict())
