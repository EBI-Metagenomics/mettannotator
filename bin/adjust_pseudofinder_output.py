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
import filecmp
import logging
import re
import shutil
import sys

logging.basicConfig(level=logging.INFO)


def main(pseudofinder_file, standard_gff, compliant_gff, outfile):
    # Check if the two GFF files are identical. If they are we don't need to do anything
    if filecmp.cmp(standard_gff, compliant_gff, shallow=False):
        logging.info("Standard and compliant GFFs are identical, no changed to Pseudofinder output needed")
        shutil.copyfile(pseudofinder_file, outfile)
    else:
        chromosome_dictionary = make_chromosome_conversion(standard_gff, compliant_gff)
        standard_genes = load_gff(standard_gff, dict())
        compliant_genes = load_gff(compliant_gff, chromosome_dictionary)
        name_translation = match_ids(standard_genes, compliant_genes)
        with open(outfile, "w") as file_out, open(pseudofinder_file, "r") as file_in:
            for line in file_in:
                if line.startswith("##sequence-region"):
                    parts = line.split(" ")
                    # replace chromosome name
                    try:
                        parts[1] = chromosome_dictionary[parts[1]]
                        line = " ".join(parts)
                    except KeyError:
                        logging.error("Could not convert sequence region {}".format(line))
                    file_out.write(line)
                elif line.startswith("#"):
                    file_out.write(line)
                else:
                    parts = line.strip().split("\t")
                    attributes_dict = dict(
                        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", parts[8])
                    )
                    for attribute_name, attribute_value in attributes_dict.items():
                        if attribute_name == "old_locus_tag":
                            locus_tags = attribute_value.split(",")
                            new_locus_tags = list()
                            for locus_tag in locus_tags:
                                if "_ign_" not in locus_tag:
                                    try:
                                        new_locus_tags.append(name_translation[locus_tag])
                                    except KeyError:
                                        logging.error(f"Locus tag {locus_tag} cannot be converted between two GFF "
                                                      f" files. Skipping")
                                else:
                                    # these are in the intergenic regions; they will be in the pseudofinder output
                                    # but not in the final genome annotation file
                                    new_locus_tags.append(locus_tag)
                            new_locus_tag_value = ",".join(new_locus_tags)
                            attributes_dict["old_locus_tag"] = new_locus_tag_value
                    new_col_9 = ""
                    for key, value in attributes_dict.items():
                        new_col_9 += f"{key}={value};"
                    parts[8] = new_col_9.rstrip(";")
                    parts[0] = chromosome_dictionary[parts[0]]
                    modified_line = "\t".join(parts)
                    file_out.write(modified_line + "\n")


def make_chromosome_conversion(standard_gff, compliant_gff):
    chromosomes_standard = get_chromosomes(standard_gff)
    chromosomes_compliant = get_chromosomes(compliant_gff)
    chromosome_dictionary = dict()
    standard_keys = list(chromosomes_standard.keys())
    compliant_keys = list(chromosomes_compliant.keys())
    for i, chr_compliant in enumerate(compliant_keys):
        chr_standard = standard_keys[i]
        if chromosomes_standard[chr_standard] == chromosomes_compliant[chr_compliant]:
            chromosome_dictionary[chr_compliant] = chr_standard
        else:
            sys.exit(f'Unable to convert chromosome names. Failed on {chr_standard} and {chr_compliant}')
    return chromosome_dictionary


def get_chromosomes(gff):
    chromosomes = dict()
    with open(gff, "r") as file_in:
        for line in file_in:
            if line.startswith("##sequence-region"):
                _, chromosome, start, end = line.strip().split(" ")
                if "|" in chromosome:
                    chromosome = chromosome.split("|")[-1]
                chromosomes[chromosome] = f'{start}-{end}'
    return chromosomes


def match_ids(standard_genes, compliant_genes):
    name_translation = {compliant_genes[key]: standard_genes[key] for key in compliant_genes if key in standard_genes}
    # log genes that are not found in both dictionaries
    unique_to_standard = set(standard_genes.keys()) - set(compliant_genes.keys())
    unique_to_compliant = set(compliant_genes.keys()) - set(standard_genes.keys())
    if len(unique_to_standard) > 0:
        logging.info("The following genes were only identified in the standard GFF:")
        for gene in unique_to_standard:
            logging.info(gene)
    if len(unique_to_compliant) > 0:
        logging.info("The following genes were only identified in the compliant GFF:")
        for gene in unique_to_compliant:
            logging.info(gene)
    return name_translation


def load_gff(gff_file, chromosome_dictionary):
    genes = dict()
    with open(gff_file, "r") as file_in:
        for line in file_in:
            if line.startswith("##FASTA"):
                return genes
            elif line.startswith("#"):
                continue
            else:
                chromosome, _, feature, start, end, _, strand, _, col9 = line.strip().split("\t")
                if "|" in chromosome:
                    chromosome = chromosome.split("|")[-1]
                if chromosome in chromosome_dictionary:
                    chromosome = chromosome_dictionary[chromosome]
                if feature == "CDS":
                    attributes_dict = dict(
                        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", col9)
                    )
                    location = "{}={}={}={}".format(chromosome, start, end, strand)
                    genes[location] = attributes_dict["locus_tag"]
    return genes


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script corrects gene IDs to account for possible discrepancies between two "
            "annotation runs (usually generated by Prokka which is run twice in mettannotator)."
        )
    )
    parser.add_argument(
        "--pseudofinder-output",
        required=True,
        help="The path to the pseudogene GFF file produced by Pseudofinder.",
    )
    parser.add_argument(
        "--standard-gff",
        required=True,
        help="The path to the GFF that contains correct gene IDs used in other annotations.",
    )
    parser.add_argument(
        "--compliant-gff",
        required=True,
        help="The path to the compliant GFF that was used for Pseudofinder.",
    )
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help="Path to the output file where the result will be saved.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.pseudofinder_output,
        args.standard_gff,
        args.compliant_gff,
        args.outfile,
    )
