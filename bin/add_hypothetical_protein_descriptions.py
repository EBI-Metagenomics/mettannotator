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


def main(ipr_types_file, ipr_file, eggnog_file, infile, outfile):
    eggnog_info = load_eggnog(eggnog_file)
    ipr_types = load_ipr_types(ipr_types_file)
    ipr_info, ipr_memberdb_only = load_ipr(ipr_file, ipr_types)
    fasta_flag = False
    with open(infile, "r") as file_in, open(outfile, "w") as file_out:
        for line in file_in:
            if not fasta_flag:
                if line.startswith("##FASTA"):
                    fasta_flag = True
                    file_out.write(line)
                elif not line.startswith("#"):
                    contig, tool, feature, start, end, blank1, strand, blank2, col9 = line.strip().split("\t")
                    if feature == "CDS":
                        attributes_dict = dict(
                            re.split(r'(?<!\\)=', item) for item in re.split(r'(?<!\\);', col9))
                        if attributes_dict["product"] == "hypothetical protein":
                            found_function, function_source = get_function(attributes_dict["ID"], attributes_dict,
                                                                           eggnog_info, ipr_info, ipr_memberdb_only)
                            found_function = escape_reserved_characters(found_function)
                            if "domain" in found_function.lower():
                                found_function = reformat_domain(found_function)
                            attributes_dict["product"] = found_function
                            attributes_dict = insert_product_source(attributes_dict, function_source)
                        else:
                            attributes_dict = insert_product_source(attributes_dict, "Prokka")
                        col9_updated = update_col9(attributes_dict)
                        file_out.write("\t".join([contig, tool, feature, start, end, blank1, strand, blank2,
                                                 col9_updated]) + "\n")
                    else:
                        file_out.write(line)
                else:
                    file_out.write(line)
            else:
                file_out.write(line)


def update_col9(attributes_dict):
    return ';'.join([f'{key}={value}' for key, value in attributes_dict.items()])


def insert_product_source(my_dict, source):
    keys_list = list(my_dict.keys())
    product_index = keys_list.index('product')
    return {k: my_dict[k] for k in keys_list[:product_index + 1]} | {"product_source": source} | \
           {k: my_dict[k] for k in keys_list[product_index + 1:]}


def get_function(acc, attributes_dict, eggnog_annot, ipr_info, ipr_memberdb_only):
    '''
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
    '''

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
        func_description, source = get_description_and_source(ipr_memberdb_only[acc], "no_type")
        if func_description:
            return func_description, source
    if acc in ipr_info and "Homologous_superfamily" in ipr_info[acc]:
        func_description, source = get_superfamily_info(ipr_info[acc]["Homologous_superfamily"])
        if func_description:
            return func_description, source
    if acc in ipr_memberdb_only and \
            any(db in {"SUPERFAMILY", "Gene3D"} for db in ipr_memberdb_only[acc]["no_type"].keys()):
        selected_keys = ['SUPERFAMILY', 'Gene3D']
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
        func_description, source = pull_out_description(my_dict[ipr_type]["NCBIfam"], "sig_desc", "ipr_desc")
        source_description = format_source("NCBIfam", source)
        return func_description, source_description
    elif any(key.lower() not in {'superfamily', 'gene3d'} for key in my_dict[ipr_type]):
        highest_match = get_best_match(my_dict[ipr_type])
        db = next(iter(highest_match))
        func_description, source = pull_out_description(highest_match[db], "ipr_desc", "sig_desc")
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
    highest_match = dict()
    for db in ipr_dict:
        if not db.lower() in ["superfamily", "gene3d"] and ipr_dict[db]["match"] > best_fraction and \
                (ipr_dict[db]["sig_desc"] != "-" or ipr_dict[db]["ipr_desc"] != "-"):
            best_fraction = ipr_dict[db]["match"]
            highest_match = dict()
            highest_match[db] = ipr_dict[db]
    if "Pfam" not in highest_match and best_fraction > 0.30:
        if "Pfam" in ipr_dict and best_fraction - ipr_dict["Pfam"]["match"] < 0.10:
            return {"Pfam": ipr_dict["Pfam"]}
    return highest_match


def find_higher_match(my_dict):
    keys = list(my_dict.keys())
    if len(keys) != 2:
        raise ValueError("The superfamily dictionary should have exactly 2 keys.")

    key1, key2 = keys
    match1 = my_dict[key1]['match']
    match2 = my_dict[key2]['match']

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
                if not cols[7] == "-" and "of unknown function" not in cols[7] and \
                        "non supervised orthologous group" not in cols[7] and "Psort location" not in cols[7] and \
                        not cols[7].lower() == "domain, protein":
                    function = cols[7]
                    # trim function from the left if it doesn't start with a letter or a digit
                    for i, char in enumerate(function):
                        if char.isalnum():
                            function = function[i:]
                            break
                    eggnog_info[cols[0]] = function
    return eggnog_info


def load_ipr_types(ipr_types_file):
    ipr_types = dict()
    with open(ipr_types_file, "r") as file_in:
        for line in file_in:
            if line.startswith("IPR"):
                acc, type, _ = line.strip().split("\t")
                ipr_types[acc] = type
    return ipr_types


def load_ipr(file, ipr_types):
    ipr_info = dict()  # hit is assigned an interpro id
    ipr_memberdb_only = dict()  # hit only exists in a member database

    with open(file, "r") as file_in:
        for line in file_in:
            cols = line.strip().split("\t")
            acc, len, db, sig_description, start, end, evalue, ipr_acc, ipr_description = \
            cols[0], cols[2], cols[3], cols[5], cols[6], cols[7], cols[8], cols[11], cols[12]
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
            if sig_description == "-" and ipr_description == "-":
                continue
            perc_match = (int(end) - int(start))/int(len)
            if perc_match < 0.10:
                continue
            if not ipr_acc == "-":
                try:
                    ipr_type = ipr_types[ipr_acc]
                except:
                    continue  # entry is no longer in InterPro
                if ipr_type not in ["Domain", "Family", "Homologous_superfamily"]:
                    continue
                ipr_info = save_to_dict(ipr_info, acc, db, perc_match, ipr_description, sig_description, ipr_type)
            else:
                ipr_memberdb_only = save_to_dict(ipr_memberdb_only, acc, db, perc_match, ipr_description,
                                                 sig_description, "no_type")
    return ipr_info, ipr_memberdb_only


def save_to_dict(res_dict, acc, db, perc_match, ipr_description, sig_description, ipr_type):
    if acc in res_dict and ipr_type in res_dict[acc] and db in res_dict[acc][ipr_type]:
        if perc_match > res_dict[acc][ipr_type][db]["match"]:
            res_dict[acc][ipr_type][db]["match"] = perc_match
            res_dict[acc][ipr_type][db]["ipr_desc"] = ipr_description
            res_dict[acc][ipr_type][db]["sig_desc"] = sig_description
    else:
        res_dict.setdefault(acc, dict())
        res_dict[acc].setdefault(ipr_type, dict())
        res_dict[acc][ipr_type].setdefault(db, dict())
        res_dict[acc][ipr_type][db]["match"] = perc_match
        res_dict[acc][ipr_type][db]["ipr_desc"] = ipr_description
        res_dict[acc][ipr_type][db]["sig_desc"] = sig_description
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
        if text[i] == ',' and text[i - 1].isdigit() and text[i + 1].isdigit():
            return True
    return False


def replace_commas(input_string):
    if ',' not in input_string:
        # If there are no commas, do nothing
        return input_string
    result = ""
    i = 0
    while i < len(input_string):
        if input_string[i] == ',':
            if is_comma_surrounded_by_digits(input_string[i - 1:i + 2]):
                result += '%2C'
            else:
                result += '/'
            i += 1  # Skip the next character as it's already processed
        else:
            result += input_string[i]
            i += 1
    return result


def reformat_domain(string):
    substrings_to_exclude = ["domain-containing", "domain contain", "domain protein",
                             "domain-related", "domain related", "domain superfamily", "domain family"]
    if all(substring not in string.lower() for substring in substrings_to_exclude):
        return string.replace("domain", "domain-containing protein")
    else:
        return string


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
        help=(
            "Path to the output file where the result will be saved."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.ipr_entries,
        args.ipr_output,
        args.eggnog_output,
        args.infile,
        args.outfile
    )
