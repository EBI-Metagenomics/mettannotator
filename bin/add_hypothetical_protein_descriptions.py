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


def main(ipr_types_file, ipr_file, hierarchy_file, eggnog_file, infile, outfile):
    eggnog_info = load_eggnog(eggnog_file)
    levels = load_hierarchy(hierarchy_file)
    ipr_types = load_ipr_types(ipr_types_file)
    ipr_info, ipr_memberdb_only, ipr_leveled_info = load_ipr(
        ipr_file, ipr_types, levels
    )

    fasta_flag = False
    with open(infile, "r") as file_in, open(outfile, "w") as file_out:
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
                        if attributes_dict["product"] == "hypothetical protein":
                            found_function, function_source = get_function(
                                attributes_dict["ID"],
                                attributes_dict,
                                eggnog_info,
                                ipr_info,
                                ipr_memberdb_only,
                            )

                            if not function_source == "UniFIRE":
                                old_func = found_function
                                found_function = clean_up_function(found_function, function_source)
                                if not old_func == found_function:
                                    print(old_func, "\t", function_source)
                                    print(found_function, "\n")
                            found_function = escape_reserved_characters(found_function)
                            attributes_dict["product"] = found_function
                            attributes_dict = insert_product_source(
                                attributes_dict, function_source
                            )
                        else:
                            attributes_dict = insert_product_source(
                                attributes_dict, "Prokka"
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


def clean_up_function(found_function, function_source):
    if "domain" in found_function.lower():
        found_function = reformat_domain(found_function)
    if found_function.lower().endswith("family") and "protein" not in found_function.lower():
        found_function = found_function + " protein"
    return found_function


def update_col9(attributes_dict):
    return ";".join([f"{key}={value}" for key, value in attributes_dict.items()])


def count_initial_dashes(s):
    return len(s) - len(s.lstrip("-"))


def load_hierarchy(parent_child_file):
    counter = dict()
    levels = dict()
    with open(parent_child_file, "r") as file_in:
        for line in file_in:
            depth = count_initial_dashes(line)
            term = line.lstrip("-").strip().split("::")[0]
            if term in levels:
                if depth > levels[term]:
                    levels[term] = depth
            else:
                levels[term] = depth
            counter[term] = counter.get(term, 0) + 1
    return levels


def insert_product_source(my_dict, source):
    keys_list = list(my_dict.keys())
    product_index = keys_list.index("product")
    return (
        {k: my_dict[k] for k in keys_list[: product_index + 1]}
        | {"product_source": source}
        | {k: my_dict[k] for k in keys_list[product_index + 1 :]}
    )


def get_function(acc, attributes_dict, eggnog_annot, ipr_info, ipr_memberdb_only):
    """
    Identify function by carrying it over from a db match. The following priority is used:
    Priority 1: UniFIRE protein recommended full name
    Priority 2: InterPro Family (NCBIfam)
    Priority 3: InterPro Family (except NCBI fam; Pfam prioritised)
    Priority 4: InterPro Domain (NCBIfam)
    Priority 5: InterPro Domain (except NCBIfam; Pfam prioritised)
    Priority 6: InterPro member db (NCBIfam)
    Priority 7: InterPro member db (except SUPERFAMILY, Gene3D, NCBIfam; Pfam prioritised)
    Priority 8: InterPro (SUPERFAMILY or Gene3D)
    Priority 9: InterPro member db (SUPERFAMILY and Gene3D)
    Priority 10: eggNOG

    :param acc: query protein accession
    :param attributes_dict: column 9 from GFF with annotation for query protein
    :param eggnog_annot: saved eggNOG annotations
    :param ipr_info: InterPro annotation
    :param ipr_memberdb_only: annotations that don't have IPR accessions but are from InterPro member databases
    :return: function, source db
    """

    if "uf_prot_rec_fullname" in attributes_dict:
        return attributes_dict["uf_prot_rec_fullname"], "UniFIRE"
    if acc in ipr_info and "Family" in ipr_info[acc]:
        func_description, source = get_description_and_source(ipr_info[acc], "Family")
        if func_description:
            return func_description, source
    if acc in ipr_info and "Domain" in ipr_info[acc]:
        func_description, source = get_description_and_source(ipr_info[acc], "Domain")
        if func_description:
            return func_description, source
    if acc in ipr_memberdb_only:
        func_description, source = get_description_and_source(
            ipr_memberdb_only[acc], "no_type"
        )
        if func_description:
            return func_description, source
    if acc in ipr_info and "Homologous_superfamily" in ipr_info[acc]:
        func_description, source = get_superfamily_info(
            ipr_info[acc]["Homologous_superfamily"]
        )
        if func_description:
            return func_description, source
    if acc in ipr_memberdb_only and any(
        db in {"SUPERFAMILY", "Gene3D"}
        for db in ipr_memberdb_only[acc]["no_type"].keys()
    ):
        selected_keys = ["SUPERFAMILY", "Gene3D"]
        subset_dict = dict()
        for key in selected_keys:
            if key in ipr_memberdb_only[acc]["no_type"]:
                subset_dict[key] = ipr_memberdb_only[acc]["no_type"][key]
        func_description, source = get_superfamily_info(subset_dict)
        if func_description:
            return func_description, source
    if acc in eggnog_annot:
        return eggnog_annot[acc], "eggNOG"
    return "hypothetical protein", "Prokka"


def get_description_and_source(my_dict, ipr_type):
    if "NCBIfam" in my_dict[ipr_type]:
        func_description, source = pull_out_description(
            my_dict[ipr_type]["NCBIfam"], "sig_desc", "ipr_desc"
        )
        source_description = format_source("NCBIfam", source)
        return func_description, source_description
    elif any(key.lower() not in {"superfamily", "gene3d"} for key in my_dict[ipr_type]):
        highest_match = get_best_match(my_dict[ipr_type])
        db = next(iter(highest_match))
        func_description, source = pull_out_description(
            highest_match[db], "ipr_desc", "sig_desc"
        )
        source_description = format_source(db, source)
        return func_description, source_description
    else:
        return "", ""


def get_superfamily_info(my_dict):
    if len(my_dict.keys()) == 1:
        db = next(iter(my_dict.keys()))
        description, source = pull_out_description(my_dict[db], "ipr_desc", "sig_desc")
        source_description = format_source(db, source)
    else:
        db = find_higher_match(my_dict)
        description, source = pull_out_description(my_dict[db], "ipr_desc", "sig_desc")
        source_description = format_source(db, source)
    return description, source_description


def format_source(db, source):
    if source == "ipr_desc":
        return "InterPro({})".format(db)
    else:
        return db


def get_best_match(ipr_dict):
    best_fraction = 0
    best_level = 0
    highest_match = dict()
    for db in ipr_dict:
        if ipr_dict[db]["level"] is None:
            ipr_dict[db]["level"] = 0
        if not db.lower() in ["superfamily", "gene3d"] and (
            ipr_dict[db]["sig_desc"] != "-" or ipr_dict[db]["ipr_desc"] != "-"
        ):
            if ipr_dict[db]["level"] > best_level or (
                ipr_dict[db]["level"] == best_level
                and ipr_dict[db]["match"] > best_fraction
            ):
                best_fraction = ipr_dict[db]["match"]
                best_level = ipr_dict[db]["level"]
                highest_match = dict()
                highest_match[db] = ipr_dict[db]
    if "Pfam" not in highest_match and best_fraction > 0.30:
        if (
            "Pfam" in ipr_dict
            and best_fraction - ipr_dict["Pfam"]["match"] < 0.10
            and ipr_dict["Pfam"]["level"] >= best_level
        ):
            return {"Pfam": ipr_dict["Pfam"]}
    return highest_match


def find_higher_match(my_dict):
    keys = list(my_dict.keys())
    if len(keys) != 2:
        raise ValueError("The superfamily dictionary should have exactly 2 keys.")

    key1, key2 = keys
    match1 = my_dict[key1]["match"]
    match2 = my_dict[key2]["match"]

    if match1 > match2:
        return key1
    elif match2 > match1:
        return key2
    else:
        return "Gene3D"  # Return Gene3d if match values are equal


def pull_out_description(my_dict, first_priority, second_priority):
    if not my_dict[first_priority] == "-":
        return my_dict[first_priority], first_priority
    else:
        return my_dict[second_priority], second_priority


def load_eggnog(file):
    eggnog_info = dict()
    with open(file, "r") as file_in:
        for line in file_in:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                try:
                    evalue = float(cols[2])
                except:
                    continue
                if evalue > 1e-10:
                    continue
                if (
                    not cols[7] == "-"
                    and "of unknown function" not in cols[7]
                    and "non supervised orthologous group" not in cols[7]
                    and "Psort location" not in cols[7]
                    and not cols[7].lower() == "domain, protein"
                    and "may contain a frame shift" not in cols[7].lower()
                ):
                    function = cols[7]
                    # trim function from the left if it doesn't start with a letter or a digit
                    function = clean_up_eggnog_function(function)
                    for i, char in enumerate(function):
                        if char.isalnum():
                            function = function[i:]
                            break
                    eggnog_info[cols[0]] = function
    return eggnog_info


def clean_up_eggnog_function(func_description):
    for i, char in enumerate(func_description):
        if char.isalnum():
            func_description = func_description[i:]
            break
    if func_description.lower().startswith("belongs to the"):
        words = func_description.split()
        func_description = " ".join(words[3:])
    func_description = func_description.replace("'phage'", "phage").replace('"phage"', "phage")
    if func_description.endswith("ase activity") and len(func_description.split()) < 5:
        func_description = func_description.replace("ase activity", "ase-like protein")
    return func_description


def load_ipr_types(ipr_types_file):
    ipr_types = dict()
    with open(ipr_types_file, "r") as file_in:
        for line in file_in:
            if line.startswith("IPR"):
                acc, type, _ = line.strip().split("\t")
                ipr_types[acc] = type
    return ipr_types


def load_ipr(file, ipr_types, ipr_levels):
    ipr_info = dict()  # hit is assigned an interpro id
    ipr_leveled_info = dict()
    ipr_memberdb_only = dict()  # hit only exists in a member database

    with open(file, "r") as file_in:
        for line in file_in:
            cols = line.strip().split("\t")
            (
                acc,
                len,
                db,
                sig_description,
                start,
                end,
                evalue,
                ipr_acc,
                ipr_description,
            ) = (
                cols[0],
                cols[2],
                cols[3],
                cols[5],
                cols[6],
                cols[7],
                cols[8],
                cols[11],
                cols[12],
            )
            if evalue == "-":
                evalue = 1
            else:
                evalue = float(evalue)
            if evalue > 1e-10:
                continue
            if db in ["ProSiteProfiles", "Coils", "MobiDBLite", "PRINTS"]:
                continue
            if sig_description.lower == "uncharacterized":
                sig_description = "-"
            if ipr_description.lower == "uncharacterized":
                ipr_description = "-"
            if db == "PANTHER":
                sig_description = clean_panther(sig_description)
            if sig_description == "-" and ipr_description == "-":
                continue
            perc_match = (int(end) - int(start)) / int(len)
            if perc_match < 0.10:
                continue
            if not ipr_acc == "-":
                if ipr_acc not in ipr_levels:
                    level = 0
                else:
                    level = ipr_levels[ipr_acc]
            if not ipr_acc == "-":
                try:
                    ipr_type = ipr_types[ipr_acc]
                except ValueError:
                    continue  # entry is no longer in InterPro
                if ipr_type not in ["Domain", "Family", "Homologous_superfamily"]:
                    continue
                ipr_info = save_to_dict(
                    ipr_info,
                    acc,
                    db,
                    perc_match,
                    ipr_description,
                    sig_description,
                    ipr_type,
                    level,
                )
            else:
                ipr_memberdb_only = save_to_dict(
                    ipr_memberdb_only,
                    acc,
                    db,
                    perc_match,
                    ipr_description,
                    sig_description,
                    "no_type",
                    None,
                )
    return ipr_info, ipr_memberdb_only, ipr_leveled_info


def save_to_dict(
    res_dict, acc, db, perc_match, ipr_description, sig_description, ipr_type, level
):
    entry = res_dict.setdefault(acc, {}).setdefault(ipr_type, {}).setdefault(db, {})

    if "level" in entry and type(level) == int:
        # For Family and Domain entries, prioritise replacement with a lower level term rather than a better percent
        # match. Lower level has a higher numerical value (level 2 is lower than 0). Only replace with a better
        # percent match of the level is not lower.
        if (
            ipr_type in ["Family", "Domain"] and level > entry["level"]
        ) or (perc_match > entry["match"] and level == entry["level"]):
            entry.update(
                {
                    "match": perc_match,
                    "ipr_desc": ipr_description,
                    "sig_desc": sig_description,
                    "level": level,
                }
            )
    else:
        entry.update(
            {
                "match": perc_match,
                "ipr_desc": ipr_description,
                "sig_desc": sig_description,
                "level": level,
            }
        )

    return res_dict


def escape_reserved_characters(string):
    string = replace_commas(string)
    reserved_characters = [";", "=", "&"]
    for ch in reserved_characters:
        if ch in string:
            string = string.replace(ch, "\{}".format(ch))
    return string


def is_comma_surrounded_by_digits(text):
    for i in range(1, len(text) - 1):
        if text[i] == "," and text[i - 1].isdigit() and text[i + 1].isdigit():
            return True
    return False


def replace_commas(input_string):
    if "," not in input_string:
        # If there are no commas, do nothing
        return input_string
    result = ""
    i = 0
    while i < len(input_string):
        if input_string[i] == ",":
            if is_comma_surrounded_by_digits(input_string[i - 1 : i + 2]):
                result += "%2C"
            else:
                result += "/"
            i += 1  # Skip the next character as it's already processed
        else:
            result += input_string[i]
            i += 1
    return result


def reformat_domain(string):
    substrings_to_exclude = [
        "domain-containing",
        "domain contain",
        "domain protein",
        "domain-related",
        "domain related",
        "domain superfamily",
        "domain family",
        "domain-",
        "protein"
    ]
    if all(substring not in string.lower() for substring in substrings_to_exclude):
        return string.replace("domain", "domain-containing protein")
    else:
        return string


def clean_panther(sig_description):
    starts_to_exclude = [
        "meiotically up-regulated gene",
    ]
    full_strings_to_exlude = [
        "family protein, putative-related",
        "putative-related",
        "protein, putative-related",
        "conserved protein",
        "domain-containing protein, putative-related",
        "unnamed product",
        "unnamed product-related"
    ]
    for start in starts_to_exclude:
        if sig_description.lower().startswith(start):
            return "-"
    for full_string in full_strings_to_exlude:
        if sig_description.lower() == full_string:
            return "-"
    return sig_description


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script Uses UniFIRE information, InterProScan and eggNOG outputs to, where possible, change "
            "product description from 'hypothetical protein' to a putative function based on homology to db entrees."
        )
    )
    parser.add_argument(
        "--ipr-entries",
        required=True,
        help="The path to the entries.list file from InterPro.",
    )
    parser.add_argument(
        "--ipr-hierarchy",
        required=True,
        help="The path to the ParentChildTreeFile.txt file from InterPro.",
    )
    parser.add_argument(
        "--ipr-output",
        required=True,
        help="The path to the TSV file produced by InterProScan.",
    )
    parser.add_argument(
        "--eggnog-output",
        required=True,
        help="The path to the TSV annotations file produced by emapper.",
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
        args.ipr_output,
        args.ipr_hierarchy,
        args.eggnog_output,
        args.infile,
        args.outfile,
    )
