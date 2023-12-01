#!/usr/bin/env python3

import argparse
import logging
import re
import sys

logging.basicConfig(level=logging.INFO)


def main(arba, unirule, pirsr, gff, outfile):
    arba_dict = load_unirule_arba(arba)
    unirule_dict = load_unirule_arba(unirule)
    pirsr_dict = load_pirsr(pirsr)
    combine_and_print(arba_dict, unirule_dict, pirsr_dict, gff, outfile)
    """
    for record in unirule_dict:
        print("\n\n{}".format(record))
        print("UNIRULE:", unirule_dict[record])
        if record in arba_dict:
            print("ARBA:", arba_dict[record])
        if record in pirsr_dict:
            print("PIRSR:", pirsr_dict[record])
    for record in arba_dict:
        if record not in unirule_dict:
            print("\n\n{}".format(record))
            print("ARBA:", arba_dict[record])
            if record in pirsr_dict:
                print("PIRSR:", pirsr_dict[record])
    for record in pirsr_dict:
        if record not in unirule_dict and record not in arba_dict:
            print("\n\n{}".format(record))
            print("PIRSR:", pirsr_dict[record])
    """

# unifire_ec_number=
# unifire_ontology_term=
# unifire_chebi=
# unifire_keywords=
# unifire_recommended_name=
# unifire_alternative_name=

def combine_and_print(arba_dict, unirule_dict, pirsr_dict, gff, outfile):
    fasta_flag = False
    list_of_dicts = [arba_dict, unirule_dict, pirsr_dict]
    combined_keys = list()
    with open(outfile, "w") as file_out, open(gff, "r") as file_in:
        for line in file_in:
            if fasta_flag is True:
                file_out.write(line)
            else:
                if line.startswith("#"):
                    file_out.write(line)
                elif line.startswith(">"):
                    file_out.write(line)
                    fasta_flag = True
                else:
                    # temp line below
                    file_out.write(line)
                    contig, tool, feature, start, end, blank1, strand, blank2, col9 = line.strip().split("\t")
                    if feature == "CDS":
                        id = get_id(col9)
                        #print("\n{}".format(id))
                        combined_dict = dict()
                        for db_dict in list_of_dicts:
                            if id in db_dict:
                                for key, value in db_dict[id].items():
                                    if key in combined_dict:
                                        combined_dict[key].extend(value)
                                    else:
                                        combined_dict[key] = value.copy()
                        if len(combined_dict) > 0:
                            combined_dict = condense_dict(combined_dict)
                            #combined_dict = escape_reserved_characters(combined_dict)
                        if combined_dict and "xref.GO" in combined_dict:
                            #print("YES", combined_dict)
                            combined_dict["xref.GO"] = ",".join(combined_dict["xref.GO"])
                            #print(combined_dict)
                        for key, value in combined_dict.items():
                            #print(key, value)
                            combined_keys.append(key)
    combined_keys = set(combined_keys)
    #print(combined_keys)


def escape_reserved_characters(combined_dict):
    reserved_characters = [";", "=", "&", ","]
    for key, list_of_values in combined_dict.items():
        for value in list_of_values:
            for ch in reserved_characters:
                if ch in value:
                    #print(key, value)
                    pass


def condense_dict(combined_dict):
    for key, list_to_condense in combined_dict.items():
        #print("Before: {}{}".format(key, list_to_condense))
        list_to_condense = list(set(list_to_condense))
        list_to_condense = collapse_keywords(list_to_condense)
        combined_dict[key] = list_to_condense
        #print("After: {}{}".format(key, list_to_condense))
    return combined_dict


def get_id(col9):
    id_pattern = r"ID=(.*?);"
    id_match = re.search(id_pattern, col9).group(1)
    return id_match


def load_pirsr(pirsr):
    results_dict = dict()
    with open(pirsr, "r") as file_in:
        for line in file_in:
            if not line.startswith("Evidence"):
                if any(keyword in line for keyword in ["keyword", "comment.cofactor"]):
                    evidence, protein_id, annot_type, value = line.strip().split("\t")
                    if annot_type == "comment.cofactor":
                        elements = value.split(";")
                        for element in elements:
                            if element.lstrip().startswith("Xref"):
                                value = element.split("=")[1]
                                annot_type = "chebi"
                            elif element.lstrip().startswith("Name"):
                                value = element.split("=")[1]
                                annot_type = "pirsr_name"
                            results_dict.setdefault(protein_id, dict()).setdefault(annot_type, list()).append(value)
                    else:
                        results_dict.setdefault(protein_id, dict()).setdefault(annot_type, list()).append(value)
    for record in results_dict:
        if "keyword" in results_dict[record]:
            if len(results_dict[record]["keyword"]) > 1:
                results_dict[record]["keyword"] = collapse_keywords(list(results_dict[record]["keyword"]))
        if "chebi" in results_dict[record]:
            if len(results_dict[record]["chebi"]) > 1:
                results_dict[record]["chebi"] = list(set(results_dict[record]["chebi"]))
    return results_dict


def load_unirule_arba(file):
    results_dict = dict()
    with open(file, "r") as file_in:
        for line in file_in:
            if not line.startswith("Evidence"):
                parts = line.strip().split("\t")
                if len(parts) == 6:
                    #evidence, protein_id, annot_type, value, start, end = line.strip().split("\t")
                    pass
                else:
                    evidence, protein_id, annot_type, value = line.strip().split("\t")
                if not annot_type.startswith("comment"):
                    results_dict.setdefault(protein_id, dict()).setdefault(annot_type, list()).append(value)
    for record in results_dict:
        if "keyword" in results_dict[record]:
            if len(results_dict[record]["keyword"]) > 1:
                results_dict[record]["keyword"] = collapse_keywords(list(results_dict[record]["keyword"]))
    return results_dict


def collapse_keywords(keyword_list):
    #print("\nKeywords before", keyword_list)
    collapsed_list = list()

    for keyword in keyword_list:
        redundant = False
        for other_keyword in keyword_list:
            if keyword != other_keyword and keyword.lower() in other_keyword.lower():
                redundant = True
        if not redundant:
            collapsed_list.append(keyword)
    collapsed_list = list(set(collapsed_list))
    #print("Keywords after", collapsed_list)
    return collapsed_list


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Placeholder"
        )
    )
    parser.add_argument(
        "-a",
        dest="arba",
        required=True,
        help="Arba.",
    )
    parser.add_argument(
        "-u",
        dest="unirule",
        required=True,
        help=(
            "Unirule."
        ),
    )
    parser.add_argument(
        "-p",
        dest="pirsr",
        required=True,
        help=(
            "Pirsr."
        ),
    )
    parser.add_argument(
        "-g",
        dest="gff",
        required=True,
        help=(
            "gff."
        ),
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help=(
            "outfile"
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.arba,
        args.unirule,
        args.pirsr,
        args.gff,
        args.outfile
    )
