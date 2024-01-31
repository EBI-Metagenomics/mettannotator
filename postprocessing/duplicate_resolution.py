#!/usr/bin/env python3

import argparse
import logging
import re

logging.basicConfig(level=logging.INFO)


def main(reference, target, outfile, species):
    # Choose stable ID prefix depending on the species
    species_prefix = get_prefix(species)
    # Load known gene names and prefixes from the reference
    alias_genename_dict = load_ref_genenames(reference, species_prefix)
    # Get duplicate gene names from the GFF we are modifying.
    # In dedupl_dict, key = base gene name (e.g. susC), value = dictionary where:
    # key = full gene name (e.g. susC_1), value = dictionary where:
    # key = "locus", value = locus name (e.g. BU_ATCC8492_00006); key = "alias", value = alias name (e.g. BACUNI_00053)
    dedupl_dict = load_duplicates(target)

    counter = 0
    stats_dict = dict()
    replacements = dict()
    for base in dedupl_dict:
        if len(dedupl_dict[base]) > 1:
            counter += 1
            print(dedupl_dict[base])
            unknown_counter = 0
            decision_dict = dict()
            for gene in dedupl_dict[base]:
                if dedupl_dict[base][gene]["alias"] in alias_genename_dict:
                    print(
                        "target",
                        dedupl_dict[base][gene]["alias"],
                        "ref",
                        alias_genename_dict[dedupl_dict[base][gene]["alias"]],
                    )
                    ref_gene_name = alias_genename_dict[dedupl_dict[base][gene]["alias"]]
                    decision_dict[ref_gene_name] = decision_dict.get(ref_gene_name, 0) + 1
                else:
                    print("target", dedupl_dict[base][gene]["alias"], "ref", "Unknown gene")
                    unknown_counter += 1
            if len(decision_dict) == 0 and unknown_counter > 1:
                stats_dict["unknowns_only"] = stats_dict.get("unknowns_only", 0) + 1
            elif any(count > 1 for count in decision_dict.values()):
                stats_dict["unable_to_decide"] = stats_dict.get("unable_to_decide", 0) + 1
                print("----------------------- > Unable to decide", decision_dict)
            else:
                print("====!====!====!====!====!=== > Can decide", decision_dict)
                if base in decision_dict and len(decision_dict) == 1:
                    print("Deciding now", dedupl_dict[base])

            # elif base in decision_dict and decision_dict[base] == 1:
            #    if len(decision_dict) == 1:
            #        stats_dict["deduplicated"] = stats_dict.get("deduplicated", 0) + 1
            #        print("========================> Should be able to decide", decision_dict)
            #    else:
            #        print("############### multiple gene names", decision_dict)
            # elif unknown_counter > 0 and len(decision_dict) == 1:
            #     ref_gene_name = list(decision_dict.keys())[0]
            #     if ref_gene_name == base:
            #         if decision_dict[ref_gene_name] == 1:
            #             stats_dict["unknowns and one matching gene"] = \
            #                 stats_dict.get("unknowns and one matching gene", 0) + 1
            #         else:
            #             stats_dict["unknowns and several matching genes"] = \
            #                 stats_dict.get("unknowns and several matching genes", 0) + 1
            #     else:
            #         stats_dict["unknowns and different genes"] = stats_dict.get("unknowns and different genes", 0) + 1
            # elif unknown_counter > 0 and len(decision_dict) > 1:
            #     stats_dict["unknowns and matching and mismatching genes"] = stats_dict.get("unknowns and matching and mismatching genes", 0) + 1
            # elif unknown_counter == 0:
            #     if ref_gene_name == base:
            #         if decision_dict[ref_gene_name] == 1:
            #             stats_dict["one matching gene"] = \
            #                 stats_dict.get("one matching gene", 0) + 1
            #         else:
            #             stats_dict["several copies match the gene name"] = \
            #                 stats_dict.get("several copies match the gene name", 0) + 1
            #     else:
            #         stats_dict["different genes"] = stats_dict.get("different genes", 0) + 1
            # else:
            #     print("------------------>UNACCOUNTED")
            #     stats_dict["unaccounted"] = \
            #         stats_dict.get("unaccounted", 0) + 1

    print("Total number of groups: {}".format(counter))
    print(stats_dict)


def load_duplicates(infile):
    dedupl_dict = dict()
    with open(infile, "r") as file_in:
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
                    except:
                        gene_name = None
                    try:
                        locus_name = re.search(locus_pattern, annot).group(1)
                    except:
                        locus_name = None
                    try:
                        alias_name = re.search(alias_pattern, annot).group(1)
                    except:
                        alias_name = None
                    print(gene_name, locus_name, alias_name)
                    if "_" in gene_name:
                        base, copy_num = gene_name.split("_")
                        try:
                            int(copy_num)
                            dedupl_dict.setdefault(base, dict()).setdefault(gene_name, dict())
                            dedupl_dict[base][gene_name] = {"locus": locus_name, "alias": alias_name}
                        except ValueError:
                            pass
    return dedupl_dict


def get_prefix(species):
    if species.lower() == "uniformis":
        return "BACUNI"
    elif species.lower() == "vulgatus":
        return "BVU"
    else:
        raise ValueError("Unknown species")


def load_ref_genenames(reference, species_prefix):
    alias_genename_dict = dict()
    with open(reference, "r") as file_in:
        for line in file_in:
            if not line.startswith("#"):
                _, _, feature, _, _, _, _, _, annot = line.strip().split("\t")
                if feature == "gene":
                    if ";Name=" in annot and species_prefix in annot:
                        id_pattern = r"ID=gene:(.*?);"
                        name_pattern = r";Name=(.*?);"
                        id = re.search(id_pattern, annot).group(1)
                        gene_name = re.search(name_pattern, annot).group(1)
                        alias_genename_dict[id] = gene_name
    return alias_genename_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script uses liftoff pipeline output to determine which of the genes that Prokka "
            "gave identical names to (for example, susC_1, susC_2, susC_3, etc) is the true one. "
            "If the script is able to determine that, it takes off the underscore. "
        )
    )
    parser.add_argument(
        "-r",
        dest="reference",
        required=True,
        help="The reference GFF files that has the 'correct' gene names.",
    )
    parser.add_argument(
        "-t",
        dest="target",
        required=True,
        help=("The GFF file generated by Prokka."),
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help=("The name of the file the result will be saved to."),
    )
    parser.add_argument(
        "-s",
        dest="species",
        required=True,
        choices=["uniformis", "vulgatus"],
        help=("The name of the species - uniformis or vulgatus."),
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
