#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright 2023 EMBL - European Bioinformatics Institute
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
import sys


def get_iprs(ipr_annot):
    iprs = {}
    with open(ipr_annot, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            protein = cols[0]
            try:
                evalue = float(cols[8])
            except:
                continue
            if evalue > 1e-10:
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
    return iprs


def get_eggnog(eggnog_annot):
    eggnogs = {}
    with open(eggnog_annot, "r") as f:
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
                try:
                    cog = cols[eggnog_fields["cog_func"]]
                    cog = list(cog)
                    if len(cog) > 1:
                        cog = ["R"]
                except Exception:
                    cog = ["NA"]
                kegg = cols[eggnog_fields["KEGG_ko"]].split(",")
                go = cols[eggnog_fields["GOs"]]
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
    # save positions of each BGC cluster to dictionary cluster_positions
    # and save the annotations to dictionary bgc_result
    with open(bgc_file, "r") as bgc_in:
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
    with open(prokka_gff, "r") as gff_in:
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
    with open(amr_file, "r") as f:
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
                    "AMRFinderPlus_gene_symbol={}".format(gene_name),
                    "AMRFinderPlus_sequence_name={}".format(seq_name),
                    "AMRFinderPlus_scope={}".format(scope),
                    "element_type={}".format(element_type),
                    "element_subtype={}".format(element_subtype),
                    "drug_class={}".format(drug_class),
                    "drug_subclass={}".format(drug_subclass),
                ]
            )
    return amr_annotations


def get_dbcan(dbcan_file):
    dbcan_annotations = dict()
    substrates = dict()
    with open(dbcan_file, "r") as f:
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
    with open(df_file, "r") as f:
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


def add_gff(
    in_gff,
    eggnog_file,
    ipr_file,
    sanntis_file,
    amr_file,
    antismash_file,
    gecco_file,
    dbcan_file,
    defense_finder_file,
):
    eggnogs = get_eggnog(eggnog_file)
    iprs = get_iprs(ipr_file)
    sanntis_bgcs = get_bgcs(sanntis_file, in_gff, tool="sanntis")
    gecco_bgcs = get_bgcs(gecco_file, in_gff, tool="gecco")
    antismash_bgcs = get_bgcs(antismash_file, in_gff, tool="antismash")
    amr_annotations = {}
    if amr_file:
        amr_annotations = get_amr(amr_file)
    dbcan_annotations = get_dbcan(dbcan_file)
    defense_finder_annotations = get_defense_finder(defense_finder_file)
    added_annot = {}
    out_gff = []
    with open(in_gff, "r") as f:
        for line in f:
            line = line.strip()
            if line[0] != "#":
                cols = line.split("\t")
                if len(cols) == 9:
                    annot = cols[8]
                    protein = annot.split(";")[0].split("=")[-1]
                    added_annot[protein] = {}
                    try:
                        eggnogs[protein]
                        pos = 0
                        for a in eggnogs[protein]:
                            pos += 1
                            if a != [""] and a != ["NA"]:
                                if pos == 1:
                                    added_annot[protein]["eggNOG"] = a
                                elif pos == 2:
                                    added_annot[protein]["COG"] = a
                                elif pos == 3:
                                    added_annot[protein]["KEGG"] = a
                                elif pos == 4:
                                    added_annot[protein]["Ontology_term"] = a
                    except Exception:
                        pass
                    try:
                        iprs[protein]
                        pos = 0
                        for a in iprs[protein]:
                            pos += 1
                            a = list(a)
                            if a != [""] and a:
                                if pos == 1:
                                    added_annot[protein]["Pfam"] = a
                                elif pos == 2:
                                    added_annot[protein]["InterPro"] = a
                    except Exception:
                        pass
                    try:
                        sanntis_bgcs[protein]
                        for key, value in sanntis_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except Exception:
                        pass
                    try:
                        gecco_bgcs[protein]
                        for key, value in gecco_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except Exception:
                        pass
                    try:
                        antismash_bgcs[protein]
                        for key, value in antismash_bgcs[protein].items():
                            added_annot[protein][key] = value
                    except Exception:
                        pass
                    try:
                        amr_annotations[protein]
                        added_annot[protein]["AMR"] = amr_annotations[protein]
                    except Exception:
                        pass
                    try:
                        dbcan_annotations[protein]
                        added_annot[protein]["dbCAN"] = dbcan_annotations[protein]
                    except Exception:
                        pass
                    try:
                        defense_finder_annotations[protein]
                        added_annot[protein]["defense_finder"] = (
                            defense_finder_annotations[protein]
                        )
                    except Exception:
                        pass
                    for a in added_annot[protein]:
                        value = added_annot[protein][a]
                        if type(value) is list:
                            value = ",".join(value)
                        if a in ["AMR", "dbCAN", "defense_finder"]:
                            cols[8] = "{};{}".format(cols[8], value)
                        else:
                            if not value == "-":
                                cols[8] = "{};{}={}".format(cols[8], a, value)
                    line = "\t".join(cols)
            out_gff.append(line)
    return out_gff


def get_rnas(ncrnas_file):
    ncrnas = {}
    counts = 0
    with open(ncrnas_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                cols = line.strip().split()
                counts += 1
                contig = cols[3]
                locus = "{}_ncRNA{}".format(contig, counts)
                product = " ".join(cols[26:])
                model = cols[2]
                strand = cols[11]
                if strand == "+":
                    start = int(cols[9])
                    end = int(cols[10])
                else:
                    start = int(cols[10])
                    end = int(cols[9])
                ncrnas.setdefault(contig, list()).append(
                    [locus, start, end, product, model, strand]
                )
                # if contig not in ncrnas:
                #    ncrnas[contig] = [[locus, start, end, product, model, strand]]
                # else:
                #    ncrnas[contig].append([locus, start, end, product, model, strand])
    return ncrnas


def load_crispr(crispr_file):
    crispr_annotations = dict()
    with open(crispr_file, "r") as f:
        for line in f:
            if not line.startswith("#"):
                contig = line.split("\t")[0]
                crispr_annotations.setdefault(contig, list())
                crispr_annotations[contig].append(line)
    return crispr_annotations


def add_ncrnas_and_crispr_to_gff(gff_outfile, ncrnas, crispr_annotations, res):
    gff_out = open(gff_outfile, "w")
    added_ncrnas = set()
    added_crisprs = set()
    for line in res:
        cols = line.strip().split("\t")
        if line[0] != "#" and len(cols) == 9:
            if cols[2] == "CDS":
                contig = cols[0]
                if contig in ncrnas:
                    for c in ncrnas[contig]:
                        locus = c[0]
                        start = str(c[1])
                        end = str(c[2])
                        product = c[3]
                        model = c[4]
                        strand = c[5]
                        if locus not in added_ncrnas:
                            added_ncrnas.add(locus)
                            annot = [
                                "ID=" + locus,
                                "inference=Rfam:14.9",
                                "locus_tag=" + locus,
                                "product=" + product,
                                "rfam=" + model,
                            ]
                            annot = ";".join(annot)
                            newLine = [
                                contig,
                                "INFERNAL:1.1.4",
                                "ncRNA",
                                start,
                                end,
                                ".",
                                strand,
                                ".",
                                annot,
                            ]
                            gff_out.write("\t".join(newLine) + "\n")
                if contig in crispr_annotations:
                    for crispr_line in crispr_annotations[contig]:
                        crispr_parts = crispr_line.strip().split("\t")
                        if (
                            "{}_{}_{}".format(
                                crispr_parts[2], crispr_parts[3], crispr_parts[4]
                            )
                            not in added_crisprs
                        ):
                            added_crisprs.add(
                                "{}_{}_{}".format(
                                    crispr_parts[2], crispr_parts[3], crispr_parts[4]
                                )
                            )
                            gff_out.write(crispr_line)
                gff_out.write("{}\n".format(line))
        #            else:
        #                gff_out.write("{}\n".format(line))
        else:
            gff_out.write("{}\n".format(line))
    gff_out.close()


if __name__ == "__main__":
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
        required=True,
    )
    parser.add_argument(
        "-e",
        dest="eggnog",
        help="eggnog annotations for the cluster repo",
        required=False,
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
    parser.add_argument("-r", dest="rfam", help="Rfam results", required=True)
    parser.add_argument("-o", dest="outfile", help="Outfile name", required=False)

    args = parser.parse_args()

    gff = args.gff_input

    extended_gff = add_gff(
        in_gff=gff,
        eggnog_file=args.eggnog,
        ipr_file=args.ips,
        sanntis_file=args.sanntis,
        amr_file=args.amr,
        antismash_file=args.antismash,
        gecco_file=args.gecco,
        dbcan_file=args.dbcan,
        defense_finder_file=args.defense_finder,
    )

    ncRNAs = get_rnas(args.rfam)
    crispr_annotations = {}
    if args.crispr:
        crispr_annotations = load_crispr(args.crispr)

    outfile = gff.split(".gff")[0] + "_annotated.gff"
    if args.outfile:
        outfile = args.outfile

    add_ncrnas_and_crispr_to_gff(outfile, ncRNAs, crispr_annotations, extended_gff)
