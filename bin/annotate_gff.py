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
import csv
import re
import sys


def main(
    gff,
    ipr_file,
    eggnog_file,
    sanntis_file,
    crispr_file,
    amr_file,
    antismash_file,
    gecco_file,
    dbcan_file,
    defense_finder_file,
    pseudofinder_file,
    rfam_file,
    trnascan_file,
    outfile,
    pseudogene_report_file,
):
    # load annotations and add them to existing CDS
    # here header contains leading GFF lines starting with "#",
    # main_gff_extended is a dictionary that contains Prokka GFF lines with added in additional annotations
    # fasta is the fasta portion of the Prokka GFF files
    header, main_gff_extended, fasta, pseudogene_report_dict = load_annotations(
        gff,
        eggnog_file,
        ipr_file,
        sanntis_file,
        amr_file,
        antismash_file,
        gecco_file,
        dbcan_file,
        defense_finder_file,
        pseudofinder_file,
    )

    ncrnas = get_ncrnas(rfam_file)
    trnas = get_trnas(trnascan_file)
    crispr_annotations = {}
    if args.crispr:
        crispr_annotations = load_crispr(crispr_file)

    write_results_to_file(
        outfile, header, main_gff_extended, fasta, ncrnas, trnas, crispr_annotations
    )

    if pseudogene_report_file:
        print_pseudogene_report(pseudogene_report_dict, pseudogene_report_file)


def write_results_to_file(
    outfile, header, main_gff_extended, fasta, ncrnas, trnas, crispr_annotations
):
    with open(outfile, "w") as file_out:
        file_out.write("\n".join(header) + "\n")
        contig_list = list(main_gff_extended.keys())
        # check if there are any contigs that don't have CDS; if so add them in
        contig_list = check_for_additional_keys(
            ncrnas, trnas, crispr_annotations, contig_list
        )
        for contig in contig_list:
            sorted_pos_list = sort_positions(
                contig, main_gff_extended, ncrnas, trnas, crispr_annotations
            )
            for pos in sorted_pos_list:
                for my_dict in (ncrnas, trnas, crispr_annotations, main_gff_extended):
                    if contig in my_dict and pos in my_dict[contig]:
                        for line in my_dict[contig][pos]:
                            if type(line) is str:
                                file_out.write(f"{line}\n")
                            else:
                                for element in line:
                                    file_out.write(element)
        for line in fasta:
            file_out.write(f"{line}\n")


def print_pseudogene_report(pseudogene_report_dict, pseudogene_report_file):
    with open(pseudogene_report_file, "w") as file_out:
        writer = csv.writer(file_out, delimiter="\t", lineterminator="\n")
        # Print header
        writer.writerow([
            'ID',
            'Pseudogene according to Bakta/Prokka',
            'Pseudogene according to Pseudofinder',
            'AntiFam hit'
        ])

        all_keys = ["gene_caller", "pseudofinder", "antifams"]
        for protein, attributes in pseudogene_report_dict.items():
            # Fill in missing attributes with False
            line = [protein] + [str(attributes.get(key, False)) for key in all_keys]
            writer.writerow(line)


def sort_positions(contig, main_gff_extended, ncrnas, trnas, crispr_annotations):
    sorted_pos_list = list()
    for my_dict in (main_gff_extended, ncrnas, trnas, crispr_annotations):
        if contig in my_dict:
            sorted_pos_list += list(my_dict[contig].keys())
    return sorted(list(set(sorted_pos_list)))


def check_for_additional_keys(ncrnas, trnas, crispr_annotations, contig_list):
    for my_dict in (ncrnas, trnas, crispr_annotations):
        dict_keys = set(my_dict.keys())
        absent_keys = dict_keys - set(contig_list)
        if absent_keys:
            contig_list = contig_list + list(absent_keys)
    return contig_list


def get_iprs(ipr_annot):
    iprs = {}
    antifams = list()
    if not ipr_annot:
        return iprs, antifams
    with open(ipr_annot) as f:
        for line in f:
            cols = line.strip().split("\t")
            protein = cols[0]
            try:
                evalue = float(cols[8])
            except ValueError:
                continue
            if evalue > 1e-10:
                continue
            if cols[3] == "AntiFam":
                antifams.append(protein)
                continue
            if protein not in iprs:
                iprs[protein] = [set(), set()]
            if cols[3] == "Pfam":
                pfam = cols[4]
                iprs[protein][0].add(pfam)
            if len(cols) > 12:
                ipr = cols[11]
                if not ipr == "-":
                    iprs[protein][1].add(ipr)
    return iprs, antifams


def get_eggnog(eggnog_annot):
    eggnogs = {}
    with open(eggnog_annot) as f:
        for line in f:
            line = line.rstrip()
            cols = line.split("\t")
            if line.startswith("#"):
                eggnog_fields = get_eggnog_fields(line)
            else:
                try:
                    evalue = float(cols[2])
                except ValueError:
                    continue
                if evalue > 1e-10:
                    continue
                protein = cols[0]
                eggnog = [cols[1]]

                cog = list(cols[eggnog_fields["cog_func"]])
                if len(cog) > 1:
                    cog = ["R"]

                kegg = cols[eggnog_fields["KEGG_ko"]].split(",")
                # Todo: I added splitting to GO, check that I don't break anything later on
                go = cols[eggnog_fields["GOs"]].split(",")
                eggnogs[protein] = [eggnog, cog, kegg, go]
    return eggnogs


def get_eggnog_fields(line):
    cols = line.strip().split("\t")
    try:
        index_of_go = cols.index("GOs")
    except ValueError:
        sys.exit("Cannot find the GO terms column.")
    if cols[8] == "KEGG_ko" and cols[15] == "CAZy":
        eggnog_fields = {"KEGG_ko": 8, "cog_func": 20, "GOs": index_of_go}
    elif cols[11] == "KEGG_ko" and cols[18] == "CAZy":
        eggnog_fields = {"KEGG_ko": 11, "cog_func": 6, "GOs": index_of_go}
    else:
        sys.exit("Cannot parse eggNOG - unexpected field order or naming")
    return eggnog_fields


def get_bgcs(bgc_file, prokka_gff, tool):
    cluster_positions = dict()
    tool_result = dict()
    bgc_annotations = dict()
    if not bgc_file:
        return bgc_annotations
    # save positions of each BGC cluster to dictionary cluster_positions
    # and save the annotations to dictionary bgc_result
    with open(bgc_file) as bgc_in:
        for line in bgc_in:
            if not line.startswith("#"):
                (
                    contig,
                    _,
                    feature,
                    start_pos,
                    end_pos,
                    _,
                    _,
                    _,
                    annotations,
                ) = line.strip().split("\t")
                if tool == "sanntis":
                    for a in annotations.split(
                        ";"
                    ):  # go through all parts of the annotation field
                        if a.startswith("nearest_MiBIG_class="):
                            class_value = a.split("=")[1]
                        elif a.startswith("nearest_MiBIG="):
                            mibig_value = a.split("=")[1]
                elif tool == "gecco":
                    for a in annotations.split(
                        ";"
                    ):  # go through all parts of the annotation field
                        if a.startswith("Type="):
                            type_value = a.split("=")[1]
                elif tool == "antismash":
                    if feature != "CDS":
                        continue
                    for a in annotations.split(
                        ";"
                    ):  # go through all parts of the annotation field
                        if a.startswith("function="):
                            type_value = a.split("=")[1]
                # save cluster positions to a dictionary where key = contig name,
                # value = list of position pairs (list of lists)
                cluster_positions.setdefault(contig, list()).append(
                    [int(start_pos), int(end_pos)]
                )
                # save BGC annotations to dictionary where key = contig, value = dictionary, where
                # key = 'start_end' of BGC, value = dictionary, where key = feature type, value = description
                if tool == "sanntis":
                    tool_result.setdefault(contig, dict()).setdefault(
                        "_".join([start_pos, end_pos]),
                        {
                            "nearest_MiBIG_class": class_value,
                            "nearest_MiBIG": mibig_value,
                        },
                    )
                elif tool == "gecco":
                    tool_result.setdefault(contig, dict()).setdefault(
                        "_".join([start_pos, end_pos]),
                        {"bgc_type": type_value},
                    )
                elif tool == "antismash":
                    tool_result.setdefault(contig, dict()).setdefault(
                        "_".join([start_pos, end_pos]),
                        {"bgc_function": type_value},
                    )
    # identify CDSs that fall into each of the clusters annotated by the BGC tool
    with open(prokka_gff) as gff_in:
        for line in gff_in:
            if not line.startswith("#"):
                matching_interval = ""
                (
                    contig,
                    _,
                    _,
                    start_pos,
                    end_pos,
                    _,
                    _,
                    _,
                    annotations,
                ) = line.strip().split("\t")
                if contig in cluster_positions:
                    for i in cluster_positions[contig]:
                        if int(start_pos) in range(i[0], i[1] + 1) and int(
                            end_pos
                        ) in range(i[0], i[1] + 1):
                            matching_interval = "_".join([str(i[0]), str(i[1])])
                            break
                # if the CDS is in an interval, save cluster's annotation to this CDS
                if matching_interval:
                    cds_id = annotations.split(";")[0].split("=")[1]
                    if tool == "sanntis":
                        bgc_annotations.setdefault(
                            cds_id,
                            {
                                "nearest_MiBIG": tool_result[contig][matching_interval][
                                    "nearest_MiBIG"
                                ],
                                "nearest_MiBIG_class": tool_result[contig][
                                    matching_interval
                                ]["nearest_MiBIG_class"],
                            },
                        )
                    elif tool == "gecco":
                        bgc_annotations.setdefault(
                            cds_id,
                            {
                                "gecco_bgc_type": tool_result[contig][
                                    matching_interval
                                ]["bgc_type"],
                            },
                        )
                    elif tool == "antismash":
                        bgc_annotations.setdefault(
                            cds_id,
                            {
                                "antismash_bgc_function": tool_result[contig][
                                    matching_interval
                                ]["bgc_function"],
                            },
                        )
            elif line.startswith("##FASTA"):
                break
    return bgc_annotations


def get_amr(amr_file):
    amr_annotations = {}
    if not amr_file:
        return amr_annotations
    with open(amr_file) as f:
        for line in f:
            if line.startswith("Protein identifier"):
                continue
            (
                protein_id,
                _,
                _,
                _,
                _,
                gene_name,
                seq_name,
                scope,
                element_type,
                element_subtype,
                drug_class,
                drug_subclass,
                _,
            ) = line.strip().split("\t", 12)
            # don't add annotations for which we don't have a protein ID (these will only be
            # available in the AMRFinderPlus TSV file)
            if protein_id == "NA":
                continue
            # check for characters that could break GFF
            if ";" in seq_name:
                seq_name = seq_name.replace(";", ",")
            if "=" in seq_name:
                seq_name = seq_name.replace("=", " ")
            amr_annotations[protein_id] = ";".join(
                [
                    f"amrfinderplus_gene_symbol={gene_name}",
                    f"amrfinderplus_sequence_name={seq_name}",
                    f"amrfinderplus_scope={scope}",
                    f"element_type={element_type}",
                    f"element_subtype={element_subtype}",
                    f"drug_class={drug_class}",
                    f"drug_subclass={drug_subclass}",
                ]
            )
    return amr_annotations


def get_dbcan(dbcan_file):
    dbcan_annotations = dict()
    substrates = dict()
    if not dbcan_file:
        return dbcan_annotations
    with open(dbcan_file) as f:
        for line in f:
            if "predicted PUL" in line:
                annot_fields = line.strip().split("\t")[8].split(";")
                for a in annot_fields:
                    if a.startswith("ID="):
                        cgc = a.split("=")[1]
                    elif a.startswith("substrate_dbcan-pul"):
                        substrate_pul = a.split("=")[1]
                    elif a.startswith("substrate_dbcan-sub"):
                        substrate_ecami = a.split("=")[1]
                substrates.setdefault(cgc, {})["substrate_ecami"] = substrate_ecami
                substrates.setdefault(cgc, {})["substrate_pul"] = substrate_pul
            elif line.startswith("#"):
                continue
            else:
                cols = line.strip().split("\t")
                prot_type = cols[2]
                annot_fields = cols[8].split(";")
                if not prot_type == "null":
                    for a in annot_fields:
                        if a.startswith("ID"):
                            acc = a.split("=")[1]
                        elif a.startswith("protein_family"):
                            prot_fam = a.split("=")[1]
                        elif a.startswith("Parent"):
                            parent = a.split("=")[1]
                    dbcan_annotations[acc] = (
                        "dbcan_prot_type={};dbcan_prot_family={};substrate_dbcan-pul={};substrate_dbcan-sub={}".format(
                            prot_type,
                            prot_fam,
                            substrates[parent]["substrate_pul"],
                            substrates[parent]["substrate_ecami"],
                        )
                    )
    return dbcan_annotations


def get_defense_finder(df_file):
    defense_finder_annotations = dict()
    type_info = dict()
    if not df_file:
        return defense_finder_annotations
    with open(df_file) as f:
        for line in f:
            if "Anti-phage system" in line:
                annot_fields = line.strip().split("\t")[8].split(";")
                for a in annot_fields:
                    if a.startswith("ID="):
                        id = a.split("=")[1]
                    elif a.startswith("type"):
                        df_type = a.split("=")[1]
                    elif a.startswith("subtype"):
                        df_subtype = a.split("=")[1]
                type_info.setdefault(id, {})["df_type"] = df_type
                type_info.setdefault(id, {})["df_subtype"] = df_subtype
            elif "DefenseFinder" in line:
                annot_fields = line.strip().split("\t")[8].split(";")
                for a in annot_fields:
                    if a.startswith("ID="):
                        id = a.split("=")[1]
                    elif a.startswith("Parent="):
                        parent = a.split("=")[1]
                defense_finder_annotations[id] = (
                    "defense_finder_type={};defense_finder_subtype={}".format(
                        type_info[parent]["df_type"], type_info[parent]["df_subtype"]
                    )
                )
    return defense_finder_annotations


def get_pseudogenes(pseudofinder_file):
    pseudogenes = dict()
    if not pseudofinder_file:
        return pseudogenes
    with open(pseudofinder_file) as file_in:
        for line in file_in:
            if not line.startswith("#"):
                col9 = line.strip().split("\t")[8]
                attributes_dict = dict(
                    re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", col9)
                )
                if "note" in attributes_dict:
                    note = attributes_dict["note"]
                else:
                    note = ""
                if "old_locus_tag" in attributes_dict:
                    tags = attributes_dict["old_locus_tag"].split(",")
                    for tag in tags:
                        if "_ign_" not in tag:
                            pseudogenes[tag] = note
    return pseudogenes


def load_annotations(
    in_gff,
    eggnog_file,
    ipr_file,
    sanntis_file,
    amr_file,
    antismash_file,
    gecco_file,
    dbcan_file,
    defense_finder_file,
    pseudofinder_file,
):
    eggnogs = get_eggnog(eggnog_file)
    iprs, antifams = get_iprs(ipr_file)
    sanntis_bgcs = get_bgcs(sanntis_file, in_gff, tool="sanntis")
    gecco_bgcs = get_bgcs(gecco_file, in_gff, tool="gecco")
    antismash_bgcs = get_bgcs(antismash_file, in_gff, tool="antismash")
    amr_annotations = get_amr(amr_file)
    dbcan_annotations = get_dbcan(dbcan_file)
    defense_finder_annotations = get_defense_finder(defense_finder_file)
    pseudogenes = get_pseudogenes(pseudofinder_file)
    pseudogene_report_dict = dict()
    added_annot = {}
    main_gff = dict()
    header = []
    fasta = []
    fasta_flag = False
    with open(in_gff) as f:
        for line in f:
            line = line.strip()
            if line[0] != "#" and not fasta_flag:
                line = line.replace("db_xref", "Dbxref")
                line = line.replace(";note=", ";Note=")
                cols = line.split("\t")
                if len(cols) == 9:
                    contig, caller, feature, start, annot = (
                        cols[0],
                        cols[1],
                        cols[2],
                        cols[3],
                        cols[8],
                    )
                    if feature != "CDS":
                        if caller == "Bakta" and feature == "region":
                            main_gff.setdefault(contig, dict()).setdefault(
                                int(start), list()
                            ).append(line)
                            continue
                        else:
                            continue
                    protein = annot.split(";")[0].split("=")[-1]
                    if protein in antifams:
                        # Don't print to the final GFF proteins that are known to not be real
                        continue
                    added_annot[protein] = {}
                    # process pseudogenes
                    if "pseudo=true" in annot.lower():
                        # fix case
                        cols[8] = annot.replace("pseudo=True", "pseudo=true")
                        # gene is already marked as a pseudogene; log it but don't add to the annotation again
                        pseudogene_report_dict.setdefault(protein, dict())
                        pseudogene_report_dict[protein]["gene_caller"] = True
                        if protein in pseudogenes:
                            pseudogene_report_dict[protein]["pseudofinder"] = True
                        else:
                            pseudogene_report_dict[protein]["pseudofinder"] = False
                    else:
                        # gene caller did not detect this protein as a pseudogene; check if pseudofinder did
                        if protein in pseudogenes:
                            pseudogene_report_dict.setdefault(protein, dict())
                            pseudogene_report_dict[protein]["gene_caller"] = False
                            pseudogene_report_dict[protein]["pseudofinder"] = True
                            added_annot[protein]["pseudo"] = "true"
                            if pseudogenes[protein]:
                                cols[8] = add_pseudogene_to_note(
                                    pseudogenes[protein], cols[8]
                                )
                    # record antifams
                    if protein in antifams:
                        pseudogene_report_dict.setdefault(protein, dict())
                        pseudogene_report_dict[protein]["antifams"] = True
                    try:
                        eggnogs[protein]
                        pos = 0
                        for a in eggnogs[protein]:
                            pos += 1
                            if a != [""] and a != ["NA"]:
                                if pos == 1:
                                    added_annot[protein]["eggNOG"] = a
                                elif pos == 2:
                                    added_annot[protein]["cog"] = a
                                elif pos == 3:
                                    added_annot[protein]["kegg"] = a
                                elif pos == 4:
                                    added_annot[protein]["Ontology_term"] = a
                    except KeyError:
                        pass
                    try:
                        iprs[protein]
                        pos = 0
                        for a in iprs[protein]:
                            pos += 1
                            a = list(a)
                            if a != [""] and a:
                                if pos == 1:
                                    added_annot[protein]["pfam"] = sorted(a)
                                elif pos == 2:
                                    added_annot[protein]["interpro"] = sorted(a)
                    except KeyError:
                        pass
                    try:
                        sanntis_bgcs[protein]
                        for key, value in sanntis_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except KeyError:
                        pass
                    try:
                        gecco_bgcs[protein]
                        for key, value in gecco_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except KeyError:
                        pass
                    try:
                        antismash_bgcs[protein]
                        for key, value in antismash_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except KeyError:
                        pass
                    try:
                        amr_annotations[protein]
                        added_annot[protein]["AMR"] = amr_annotations[protein]
                    except KeyError:
                        pass
                    try:
                        dbcan_annotations[protein]
                        added_annot[protein]["dbCAN"] = dbcan_annotations[protein]
                    except KeyError:
                        pass
                    try:
                        defense_finder_annotations[protein]
                        added_annot[protein]["defense_finder"] = (
                            defense_finder_annotations[protein]
                        )
                    except KeyError:
                        pass
                    for a in added_annot[protein]:
                        value = added_annot[protein][a]
                        if type(value) is list:
                            value = ",".join(value)
                        if a in ["AMR", "dbCAN", "defense_finder"]:
                            cols[8] = f"{cols[8]};{value}"
                        else:
                            if not value == "-":
                                cols[8] = f"{cols[8]};{a}={value}"
                    line = "\t".join(cols)
                    main_gff.setdefault(contig, dict()).setdefault(
                        int(start), list()
                    ).append(line)
            elif line.startswith("#"):
                if line == "##FASTA":
                    fasta_flag = True
                    fasta.append(line)
                else:
                    header.append(line)
            elif fasta_flag:
                fasta.append(line)
    return header, main_gff, fasta, pseudogene_report_dict


def add_pseudogene_to_note(note_text, col9):
    col9_dict = dict(
        re.split(r"(?<!\\)=", item) for item in re.split(r"(?<!\\);", col9)
    )
    if "Note" in col9_dict.keys():
        col9_dict["Note"] = col9_dict["Note"] + f", {note_text}"
        return ";".join([f"{key}={value}" for key, value in col9_dict.items()])
    else:
        # insert note after locus tag
        keys_list = list(col9_dict.keys())
        locus_tag_index = keys_list.index("locus_tag")
        new_dict = (
            {k: col9_dict[k] for k in keys_list[: locus_tag_index + 1]}
            | {"Note": note_text}
            | {k: col9_dict[k] for k in keys_list[locus_tag_index + 1 :]}
        )
        return ";".join([f"{key}={value}" for key, value in new_dict.items()])


def get_ncrnas(ncrnas_file):
    ncrnas = {}
    counts = 0
    with open(ncrnas_file) as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split()
                counts += 1
                contig = cols[3]
                locus = f"{contig}_ncRNA{counts}"
                product = cols[-1]
                model = cols[2]
                if model == "RF00005":
                    # Skip tRNAs, we add them from tRNAscan-SE
                    continue
                strand = cols[11]
                if strand == "+":
                    start = int(cols[9])
                    end = int(cols[10])
                else:
                    start = int(cols[10])
                    end = int(cols[9])
                annot = [
                    "ID=" + locus,
                    "inference=Rfam:14.9",
                    "locus_tag=" + locus,
                    "product=" + product,
                    "rfam=" + model,
                ]
                annot = ";".join(annot)
                newline = "\t".join(
                    [
                        contig,
                        "INFERNAL:1.1.4",
                        "ncRNA",
                        str(start),
                        str(end),
                        ".",
                        strand,
                        ".",
                        annot,
                    ]
                )
                ncrnas.setdefault(contig, dict()).setdefault(start, list()).append(
                    newline
                )
    return ncrnas


def get_trnas(trnas_file):
    trnas = {}
    with open(trnas_file) as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.split("\t")
                contig, feature, start = cols[0], cols[2], cols[3]
                if feature == "tRNA":
                    line = line.replace("tRNAscan-SE", "tRNAscan-SE:2.0.9")
                    trnas.setdefault(contig, dict()).setdefault(
                        int(start), list()
                    ).append(line.strip())
    return trnas


def load_crispr(crispr_file):
    crispr_annotations = dict()
    with open(crispr_file) as f:
        record = list()
        left_coord = ""
        loc_contig = ""
        previous_end = ""
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split("\t")
                contig, _, start, end = (
                    cols[0],
                    cols[2],
                    int(cols[3]),
                    int(cols[4]),
                )
                if (
                    len(record) > 0
                    and contig == loc_contig
                    and abs(start - previous_end) < 2
                ):
                    # the line is a continuation of an existing record
                    record.append(line)
                    previous_end = end
                elif len(record) == 0:
                    record.append(line)
                    left_coord = start
                    loc_contig = contig
                    previous_end = end
                else:
                    # the previous record is complete, started reading a new record
                    crispr_annotations.setdefault(contig, dict()).setdefault(
                        left_coord, list()
                    ).append(record)
                    record = list()
                    record.append(line)
                    previous_end = end
                    left_coord = start
        if len(record) > 0:
            crispr_annotations.setdefault(contig, dict()).setdefault(
                left_coord, list()
            ).append(record)
    return crispr_annotations


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add functional annotation to GFF file",
    )
    parser.add_argument(
        "-g",
        dest="gff_input",
        required=True,
        help="GFF input file",
    )
    parser.add_argument(
        "-i",
        dest="ips",
        help="InterproScan annotations results for the cluster rep",
        required=False,
    )
    parser.add_argument(
        "-e",
        dest="eggnog",
        help="eggnog annotations for the cluster repo",
        required=True,
    )
    parser.add_argument(
        "-s",
        dest="sanntis",
        help="Sanntis results for the cluster rep",
        required=False,
    )
    parser.add_argument(
        "-c",
        dest="crispr",
        help="CRISPRCasFinder results for the cluster rep (high quality GFF)",
        required=False,
    )
    parser.add_argument(
        "-a",
        dest="amr",
        help="The TSV file produced by AMRFinderPlus",
        required=False,
    )
    parser.add_argument(
        "--antismash",
        help="The GFF file produced by AntiSMASH post-processing script",
        required=False,
    )
    parser.add_argument(
        "--gecco",
        help="The GFF file produced by GECCO",
        required=False,
    )
    parser.add_argument(
        "--dbcan",
        help="The GFF file produced by dbCAN post-processing script",
        required=False,
    )
    parser.add_argument(
        "--defense-finder",
        help="The GFF file produced by Defense Finder post-processing script",
        required=False,
    )
    parser.add_argument(
        "--pseudofinder",
        help="The GFF file produced by the Pseudofinder post-processing script",
        required=False,
    )
    parser.add_argument("-r", dest="rfam", help="Rfam results", required=True)
    parser.add_argument(
        "-t", dest="trnascan", help="tRNAScan-SE results", required=True
    )
    parser.add_argument("-o", dest="outfile", help="Outfile name", required=True)
    parser.add_argument(
        "--pseudogene-report", help="Pseudogene report filename", required=False
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.gff_input,
        args.ips,
        args.eggnog,
        args.sanntis,
        args.crispr,
        args.amr,
        args.antismash,
        args.gecco,
        args.dbcan,
        args.defense_finder,
        args.pseudofinder,
        args.rfam,
        args.trnascan,
        args.outfile,
        args.pseudogene_report,
    )
