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
import logging
import re
import sys

logging.basicConfig(level=logging.INFO)


def main(reference, target, outfile, species):
    # Choose stable ID prefix depending on the species
    species_prefix = get_prefix(species)
    # Load known gene names and aliases from the reference
    alias_genename_dict = load_ref_genenames(reference, species_prefix)
    # Get duplicate gene names from the GFF we are modifying.
    # In dedupl_dict, key = base gene name (e.g. susC), value = dictionary where:
    # key = full gene name (e.g. susC_1), value = dictionary where:
    # key = "locus", value = locus name (e.g. BU_ATCC8492_00006); key = "alias", value = alias name (e.g. BACUNI_00053)
    # gene_occurrence_counter has the number of times every single gene appears in the genome we are deduplicating
    dedupl_dict, gene_occurrence_counter, alias_repeats = load_duplicates(target)
    # Todo: articifially create a file with alias repeated and gene name assigned to test deduplication
    #  with extra_copy_number > 1

    counter = 0  # track total number of gene groups to deduplicate
    stats_dict = dict()  # stats for printing
    replacements = dict()  # gene names to change
    replacements_ids = dict()  # changes to make if alias is None
    reverse = list()  # undo some planned changes (remove these from replacements)
    stats_out = open(f"{species}_stats.txt", "w")
    stats_out.write("Gene\tCopy number\tReplaced\tWhy not replaced/ Comment\n")
    for base in dedupl_dict:
        if (
            len(dedupl_dict[base]) > 1
        ):  # there are multiple gene copies, try to deduplicate
            printed_stat = ""
            printed_stat += f"{base}\t{str(len(dedupl_dict[base]))}\t"
            counter += 1
            unknown_counter = 0
            decision_dict = dict()
            replace = False
            # lookup these aliases in the reference to see if they have gene names
            # if they don't, we consider them an unknown
            for gene in dedupl_dict[base]:
                if dedupl_dict[base][gene]["alias"] in alias_genename_dict:
                    ref_gene_name = alias_genename_dict[
                        dedupl_dict[base][gene]["alias"]
                    ]
                    decision_dict.setdefault(ref_gene_name, list()).append(
                        dedupl_dict[base][gene]["alias"]
                    )
                else:
                    unknown_counter += 1
            # Go through decision dictionary and see what we can resolve.
            # Decision dictionary has the following format (2 examples here):
            # {'fabHA': ['BACUNI_03004', 'BACUNI_04476']}
            # {'der': ['BACUNI_03006'], 'feoB': ['BACUNI_04461'], 'hydF': ['BACUNI_02997']}}
            logging.debug(f"Decision dictionary: {decision_dict}")
            if len(decision_dict) == 0:
                # unable to resolve duplicates, leave records as they are
                stats_dict["unknowns_only"] = stats_dict.get("unknowns_only", 0) + 1
                printed_stat += "No\tNo known genes in reference\n"
            elif any(len(alias_list) > 1 for alias_list in decision_dict.values()):
                # unable to decide because a gene name is repeated in reference as well
                stats_dict["unable_to_decide"] = (
                    stats_dict.get("unable_to_decide", 0) + 1
                )
                printed_stat += "No\tGene name occurs in reference multiple times or alias repeated in target\n"
            elif base in decision_dict and len(decision_dict) == 1:
                # there is one clear "real" gene
                alias = decision_dict[base][0]
                if alias_repeats[alias] > 1:
                    # The same alias is assigned to multiple genes, we cannot resolve this
                    stats_dict["unable_to_decide"] = (
                        stats_dict.get("unable_to_decide", 0) + 1
                    )
                    printed_stat += "No\tThe same alias is assigned to multiple genes\n"
                else:
                    if alias in replacements:
                        sys.exit(
                            f"Error: something went wrong, alias {alias} is already in replacements"
                        )
                    replacements[alias] = base
                    stats_dict["replaced"] = stats_dict.get("replaced", 0) + 1
                    replace = True
                    printed_stat += "Yes\t\n"
            else:
                logging.debug(f"Case is not clear {dedupl_dict[base]} {decision_dict}")
                unique = check_value_uniqueness(
                    decision_dict
                )  # check that the alias doesn't repeat
                occurrence_flag = False  # check if a gene we are replacing with is already in the genome
                for gene in decision_dict:
                    if not gene == base:
                        if gene in gene_occurrence_counter:
                            occurrence_flag = True
                if not unique:
                    # aliases repeat and we can't resolve duplicates
                    stats_dict["unable_to_decide"] = (
                        stats_dict.get("unable_to_decide", 0) + 1
                    )
                    printed_stat += "No\tAt least one alias is repeated\n"
                elif occurrence_flag:
                    # the gene is already in the genome and we will create a new duplicate if we use it
                    stats_dict["unable_to_decide"] = (
                        stats_dict.get("unable_to_decide", 0) + 1
                    )
                    printed_stat += (
                        "No\tReplacement already occurs elsewhere in the genome\n"
                    )
                    logging.debug(
                        "Replacement gene already occurs in the genome, can't replace"
                    )
                else:
                    already_present_in_replacements, reverse = check_gene_presence(
                        decision_dict, replacements, reverse
                    )

                    if already_present_in_replacements:
                        stats_dict["unable_to_decide"] = (
                            stats_dict.get("unable_to_decide", 0) + 1
                        )
                        printed_stat += "No\tWe already used the replacement gene in a previous duplicate group\n"
                        logging.debug(
                            "Replacement gene already occurs in the replacement list, can't replace"
                        )
                        continue

                    if len(dedupl_dict[base]) == (
                        len(decision_dict) + unknown_counter
                    ) and all(len(value) == 1 for value in decision_dict.values()):
                        # two dictionaries are the same length and new gene names don't occur in the genome
                        # we can replace every gene
                        for gene_name, alias_list in decision_dict.items():
                            alias = alias_list[0]
                            if alias not in replacements:
                                replacements[alias] = gene_name
                            else:
                                (f"Alias {alias} is already in replacements")
                                sys.exit()

                        logging.debug(f"Replaced one for one {decision_dict}")
                        stats_dict["replaced"] = stats_dict.get("replaced", 0) + 1
                        replace = True
                        printed_stat += "Yes\t\n"
                    else:
                        logging.debug("length is different")
                        if len(dedupl_dict[base]) - len(decision_dict) > 1:
                            stats_dict["unable_to_decide"] = (
                                stats_dict.get("unable_to_decide", 0) + 1
                            )
                            printed_stat += "No\tSource and replacement dicts have different lengths\n"
                        else:
                            replacements, reverse, stats_dict, printed_stat = (
                                resolve_duplicate(
                                    dedupl_dict[base],
                                    decision_dict,
                                    replacements,
                                    reverse,
                                    stats_dict,
                                    base,
                                    printed_stat,
                                )
                            )
                            if "Yes" in printed_stat:
                                replace = True
        if replace:
            replacements, printed_stat, replacements_ids = (
                try_to_remove_more_underscores(
                    replacements,
                    base,
                    printed_stat,
                    decision_dict,
                    dedupl_dict[base],
                    replacements_ids,
                )
            )
        stats_out.write(printed_stat)
    print("replacement ids", replacements_ids)
    if (
        len(set(replacements.values())) != len(replacements.values())
        or len(reverse) > 0
    ):
        sys.exit("Non-unique values in replacements")
    else:
        made_replacements = make_replacement_file(
            target, outfile, replacements, replacements_ids, species
        )
    if made_replacements != (len(replacements) + len(replacements_ids)):
        sys.exit(
            f"Made {str(made_replacements)} replacements but expected {str(len(replacements))}"
        )
    logging.info(f"Total number of groups: {counter}")
    logging.info(f"Replacements: {replacements}")
    logging.info(f"Made replacements: {made_replacements}")
    logging.info(f"Reverse {reverse}")
    stats_out.close()
    return stats_dict


def try_to_remove_more_underscores(
    replacements,
    base,
    printed_stat,
    decision_dict,
    deduplication_section,
    replacements_ids,
):
    # Check if after the replacements have been made we can remove more underscores.
    # This can only happen if replacements have been made to genes that are not the base value
    # For example, you start with dnaA_1, dnaA_2, dnaA_3 and replace with dnaA_1, dnaB, dnaC
    # You can now remove the underscore from dnaA_1 because it's no longer a duplicate.
    if base not in decision_dict:
        if len(deduplication_section.keys()) - len(decision_dict.keys()) == 1:
            # one name remains not replaced and it's not already being used in deduplication
            bacunis_already_replaced = [
                item for sublist in decision_dict.values() for item in sublist
            ]
            for key in deduplication_section.keys():
                if deduplication_section[key]["alias"] is None:
                    print("Alias is none", deduplication_section[key])
                    replacements_ids[deduplication_section[key]["locus"]] = key.split(
                        "_"
                    )[0]
                    printed_stat = printed_stat.replace(
                        "\n",
                        "Removed underscore from {} because it was the only duplicate "
                        "remaining after other replacements; no stable ID was assigned to this gene\n",
                    ).format(key)
                    return replacements, printed_stat, replacements_ids
                if deduplication_section[key]["alias"] not in bacunis_already_replaced:
                    replacements[deduplication_section[key]["alias"]] = key.split("_")[
                        0
                    ]
                    printed_stat = printed_stat.replace(
                        "\n",
                        "Removed underscore from {} because it was the only duplicate "
                        "remaining after other replacements\n",
                    ).format(key)
    return replacements, printed_stat, replacements_ids


def make_replacement_file(target, outfile, replacements, replacements_ids, species):
    seq_flag = False
    count_replacements = list()
    rep_out = open(f"{species}_replacements.txt", "w")
    with open(target) as file_in, open(outfile, "w") as file_out:
        for line in file_in:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    seq_flag = True
                file_out.write(line)
            elif seq_flag:
                file_out.write(line)
            else:
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    id = fields[8].split(";")[0].split("=")[1]
                    print(id)
                    alias_pattern = r";Alias=(.*?);"
                    try:
                        alias_name = re.search(alias_pattern, fields[8]).group(1)
                    except:  # noqa: E722
                        alias_name = None
                    gene_alias_name = alias_name
                    gene_id = id
                if fields[2] in ["gene", "CDS", "mRNA", "exon"]:
                    gene_pattern = r";gene=(.*?);"
                    try:
                        gene_name = re.search(gene_pattern, fields[8]).group(1)
                    except:  # noqa E722
                        gene_name = None
                    if fields[2] in ["CDS", "mRNA", "exon"]:
                        alias_name = gene_alias_name
                        id = gene_id
                    if alias_name and alias_name in replacements:
                        if fields[2] == "gene":
                            rep_out.write(f"{gene_name}\t{replacements[alias_name]}\n")
                        line = line.replace(gene_name, replacements[alias_name])
                        count_replacements.append(alias_name)
                    print(f"ID is {id}")
                    if id in replacements_ids:
                        print("ID IS IN IDS")
                        if fields[2] == "gene":
                            rep_out.write(f"{gene_name}\t{replacements_ids[id]}\n")
                            print(
                                "made replacement of", gene_name, replacements_ids[id]
                            )
                        line = line.replace(gene_name, replacements_ids[id])
                        count_replacements.append(id)
                    file_out.write(line)
                else:
                    file_out.write(line)
    rep_out.close()
    return len(set(count_replacements))


def resolve_duplicate(
    genes_to_resolve,
    reference_dict,
    replacements,
    reverse,
    stats_dict,
    base,
    printed_stat,
):
    logging.debug(f"=============>RESOLVING {genes_to_resolve} {reference_dict}")
    resolved = dict()
    duplicate_removed = False
    for gene in reference_dict:
        alias = reference_dict[gene][0]
        resolved[alias] = gene
        if gene == base:
            duplicate_removed = True
        for duplicate_gene, info_dict in genes_to_resolve.items():
            if genes_to_resolve[duplicate_gene]["alias"] == alias:
                del genes_to_resolve[duplicate_gene]
                break
    if len(genes_to_resolve) > 1:
        # unable to fully resolve
        logging.debug("Unable to fully resolve")
        stats_dict["unable_to_decide"] = stats_dict.get("unable_to_decide", 0) + 1
        printed_stat += "No\tResolve function can't resolve\n"
    elif len(genes_to_resolve) == 1 and not duplicate_removed:
        for duplicate_gene, info_dict in genes_to_resolve.items():
            resolved[genes_to_resolve[duplicate_gene]["alias"]] = base
            stats_dict["replaced"] = stats_dict.get("replaced", 0) + 1
        printed_stat += "Yes\t\n"
    else:
        sys.exit("Unknown case")
    logging.debug(f"============>Resolved {resolved}")
    return replacements, reverse, stats_dict, printed_stat


def check_gene_presence(decision_dict, replacements, reverse):
    result = False
    check_values = [
        value for values_list in decision_dict.values() for value in values_list
    ]
    for val in check_values:
        if val in replacements.values():
            result = True
            reverse.append(val)
    return result, reverse


def check_value_uniqueness(my_dict):
    all_values = [value for values in my_dict.values() for value in values]

    # Check for uniqueness
    if len(all_values) == len(set(all_values)):
        return True
    else:
        return False


def load_duplicates(infile):
    dedupl_dict = dict()
    gene_occurrence_counter = dict()
    alias_repeats = dict()
    with open(infile) as file_in:
        for line in file_in:
            if not line.startswith("#"):
                if line.startswith(">"):
                    break
                _, _, feature, _, _, _, _, _, annot = line.strip().split("\t")
                if feature == "gene" and "gene=" in annot:
                    gene_pattern = r";gene=(.*?);"
                    locus_pattern = r";locus_tag=(.*?);"
                    alias_pattern = r";Alias=(.*?);"
                    try:
                        gene_name = re.search(gene_pattern, annot).group(1)
                    except:  # noqa E722
                        gene_name = None
                    try:
                        locus_name = re.search(locus_pattern, annot).group(1)
                    except:  # noqa E722
                        locus_name = None
                    try:
                        alias_name = re.search(alias_pattern, annot).group(1)
                    except:  # noqa E722
                        alias_name = None
                    if alias_name:
                        alias_repeats[alias_name] = alias_repeats.get(alias_name, 0) + 1
                    if "_" in gene_name:
                        base, copy_num = gene_name.split("_")
                        gene_occurrence_counter[base] = (
                            gene_occurrence_counter.get(base, 0) + 1
                        )
                        try:
                            int(copy_num)
                            dedupl_dict.setdefault(base, dict()).setdefault(
                                gene_name, dict()
                            )
                            dedupl_dict[base][gene_name] = {
                                "locus": locus_name,
                                "alias": alias_name,
                            }
                        except ValueError:
                            pass
                    else:
                        gene_occurrence_counter[gene_name] = (
                            gene_occurrence_counter.get(gene_name, 0) + 1
                        )
    return dedupl_dict, gene_occurrence_counter, alias_repeats


def get_prefix(species):
    if species.lower() == "uniformis":
        return "BACUNI"
    elif species.lower() == "vulgatus":
        return "BVU"
    else:
        raise ValueError("Unknown species")


def load_ref_genenames(reference, species_prefix):
    alias_genename_dict = dict()
    with open(reference) as file_in:
        for line in file_in:
            if not line.startswith("#"):
                _, _, feature, _, _, _, _, _, annot = line.strip().split("\t")
                if feature == "gene":
                    if ";Name=" in annot and species_prefix in annot:
                        if species_prefix == "BACUNI":
                            id_pattern = r"ID=gene:(.*?);"
                            name_pattern = r";Name=(.*?);"
                        elif species_prefix == "BVU":
                            if "old_locus_tag" in annot and "gene=" in annot:
                                id_pattern = r"old_locus_tag=(.*?)$"
                                name_pattern = r";gene=(.*?);"
                            else:
                                continue
                        id = re.search(id_pattern, annot).group(1)
                        gene_name = re.search(name_pattern, annot).group(1)
                        alias_genename_dict[id] = gene_name
    print(alias_genename_dict)
    return alias_genename_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "This project-specific script uses liftoff pipeline output to determine which of the genes that Prokka "
            "gave identical names to (for example, susC_1, susC_2, susC_3, etc) is the true one. "
            "If the script is able to determine that, it takes off the underscore, otherwise gene names remain "
            "unchanged. The script is not intended for wider use outside of the specific project."
        )
    )
    parser.add_argument(
        "-r",
        dest="reference",
        required=True,
        help="The reference GFF file that has the 'correct' gene names.",
    )
    parser.add_argument(
        "-t",
        dest="target",
        required=True,
        help="The GFF file generated by Prokka with aliases added from reference by liftoff.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="The name of the file the result will be saved to.",
    )
    parser.add_argument(
        "-s",
        dest="species",
        required=True,
        choices=["uniformis", "vulgatus"],
        help="The name of the species - uniformis or vulgatus.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.reference,
        args.target,
        args.outfile,
        args.species,
    )
