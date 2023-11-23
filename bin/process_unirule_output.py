#!/usr/bin/env python3

import argparse
import logging
import re

logging.basicConfig(level=logging.INFO)


def main(arba, unirule, pirsr, gff, outfile):
    arba_dict = load_arba(arba)
    unirule_dict = load_unirule(unirule)
    # third db - keep comment.cofactor and keywords
    for record in unirule_dict:
        print("UNIRULE:", unirule_dict[record])
        if record in arba_dict:
            print("ARBA:", arba_dict[record])


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
    print("\nKeywords before", keyword_list)
    collapsed_list = list()

    for keyword in keyword_list:
        redundant = False
        for other_keyword in keyword_list:
            if keyword != other_keyword and keyword.lower() in other_keyword.lower():
                redundant = True
        if not redundant:
            collapsed_list.append(keyword)
    collapsed_list = list(set(collapsed_list))
    print("Keywords after", collapsed_list)
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
