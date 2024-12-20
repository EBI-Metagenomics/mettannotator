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
import re


def main(infile, outfile):
    with open(infile, "r") as file_in, open(outfile, "w") as file_out:
        fasta_flag = False
        for line in file_in:
            if line.startswith("##FASTA"):
                fasta_flag = True
                file_out.write(line)
            elif line.startswith("#"):
                file_out.write(line)
            elif fasta_flag:
                file_out.write(line)
            else:
                new_line = modify_line(line)
                file_out.write(new_line + "\n")


def modify_line(line):
    columns = line.strip().split("\t")
    if columns[2] in ["region", "tRNA"]:
        return line.strip()

    # Parse the 9th field into key-value pairs
    attributes_dict = dict(
        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", columns[8])
    )
    # Move all GO terms into Dbxref and deduplicate if needed
    attributes_dict = collect_go_terms(attributes_dict)

    # Move Pfam, Rfam and InterPro into Dbxref
    terms_to_move = ["pfam", "rfam", "interpro"]
    attributes_dict = move_to_dbxref(attributes_dict, terms_to_move)

    attributes_to_move_to_note = {
        "antismash_bgc_function": "antismash_bgc_function",
        "nearest_MiBIG": "nearest_MiBIG",
        "nearest_MiBIG_class": "nearest_MiBIG_class",
        "gecco_bgc_type": "gecco_bgc_type",
        "cog": "COG",
        "kegg": "KEGG",
        "amrfinderplus_gene_symbol": "amrfinderplus_gene_symbol",
        "amrfinderplus_scope": "amrfinderplus_scope",
        "amrfinderplus_sequence_name": "amrfinderplus_sequence_name",
        "element_type": "element_type",
        "element_subtype": "element_subtype",
        "dbcan_prot_family": "dbcan_protein_family",
        "dbcan_prot_type": "dbcan_protein_type",
        "drug_class": "drug_class",
        "drug_subclass": "drug_subclass",
        "eggNOG": "eggNOG",
        "substrate_dbcan-pul": "dbcan-PUL_substrate",
        "substrate_dbcan-sub": "dbcan-subfam_substrate",
        "uf_chebi": "UniFIRE_ChEBI",
        "uf_gene_name": "UniFIRE_gene_name",
        "uf_gene_name_synonym": "UniFIRE_gene_name_synonym",
        "uf_keyword": "UniFIRE_keyword",
        "uf_pirsr_cofactor": "UniFIRE_pirsr_cofactor",
        "uf_prot_alt_ecnumber": "UniFIRE_protein_alternative_EC_number",
        "uf_prot_alt_fullname": "UniFIRE_protein_alternative_fullname",
        "uf_prot_alt_shortname": "UniFIRE_protein_alternative_shortname",
        "uf_prot_rec_ecnumber": "UniFIRE_protein_recommended_EC_number",
        "uf_prot_rec_fullname": "UniFIRE_protein_recommended_fullname",
        "uf_prot_rec_shortname": "UniFIRE_protein_recommended_shortname",
        "defense_finder_type": "anti-phage_system_type",
        "defense_finder_subtype": "anti-phage_system_subtype",
    }
    attributes_dict = move_values_to_note(attributes_dict, attributes_to_move_to_note)

    # Reconstruct the 9th field
    columns[8] = ";".join([f"{key}={value}" for key, value in attributes_dict.items()])

    # Return the modified line
    return "\t".join(columns)


def move_values_to_note(attributes_dict, attributes_to_move_to_note):
    for term, reformatted_term in attributes_to_move_to_note.items():
        if term in attributes_dict:
            if "Note" in attributes_dict:
                attributes_dict[
                    "Note"
                ] += f",{reformatted_term}:{attributes_dict[term]}"
            else:
                attributes_dict["Note"] = f"{reformatted_term}:{attributes_dict[term]}"
            attributes_dict.pop(term)
    return attributes_dict


def move_to_dbxref(attributes_dict, terms_to_move):
    reformatted_term_names = {
        "pfam": "PFAM",
        "rfam": "RFAM",
        "interpro": "InterPro",
    }
    formatted_list_to_move = []
    for term in terms_to_move:
        if term in attributes_dict:
            reformatted_term_name = reformatted_term_names[term]
            formatted_list_to_move.extend(
                f"{reformatted_term_name}:{element}"
                for element in attributes_dict[term].split(",")
            )
            attributes_dict.pop(term)
    if len(formatted_list_to_move) > 0:
        text_to_add = ",".join(formatted_list_to_move)
        if "Dbxref" in attributes_dict:
            attributes_dict["Dbxref"] = attributes_dict["Dbxref"] + "," + text_to_add
        else:
            attributes_dict["Dbxref"] = text_to_add
    attributes_dict["Dbxref"] = attributes_dict["Dbxref"].rstrip(",").lstrip(",")
    return attributes_dict


def collect_go_terms(attributes_dict):
    # Function collects GO terms from across the 9th column into one list and deduplicates them
    go_terms = list()
    from_fields = ["Dbxref", "Ontology_term", "uf_ontology_term"]
    for field in from_fields:
        if field in attributes_dict:
            annotations = attributes_dict[field].split(",")
            for annotation in annotations[
                :
            ]:  # Using a copy of the list to avoid iteration issues
                if annotation.startswith("GO:"):
                    go_terms.append(annotation)
                    annotations.remove(annotation)
            if len(annotations) > 0:
                modified_annotations = ",".join(annotations)
                attributes_dict[field] = modified_annotations
            else:
                del attributes_dict[field]
    # deduplicate
    go_terms = list(set(go_terms))
    if "Dbxref" in attributes_dict:
        attributes_dict["Dbxref"] = attributes_dict["Dbxref"] + "," + ",".join(go_terms)
    else:
        attributes_dict["Dbxref"] = ",".join(go_terms)
    attributes_dict["Dbxref"] = attributes_dict["Dbxref"].rstrip(",")
    return attributes_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Script processes final GFF to prepare for conversion "
        "for ENA/GenBank."
    )
    parser.add_argument(
        "-i", "--infile", required=True, help="Path to the GFF file to plot"
    )
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the output file"
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.infile,
        args.outfile,
    )
