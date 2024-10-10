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

# ruff: noqa: F841, N806, N803

import argparse
import re


def gff_parser(current_line):
    (
        contig,
        seq_source,
        seq_type,
        start,
        end,
        score,
        strand,
        phase,
        attr,
    ) = current_line.split("\t")
    seq_id = attr.split(";")[0].replace("ID=", "")
    data_list = [contig, seq_type, start, end, seq_id]
    return data_list


def momo_parser(mobannot):
    ### Parsing the MAP gff file
    momofy_dict = {}
    momo_subtypes = {}
    bound_map_1 = {}
    bound_map_2 = {}
    irdr_map = {}
    irdr_dict = {}
    mges = [
        "insertion_sequence",
        "integron",
        "conjugative_transposon",
        "plasmid",
        "viral_sequence",
        "prophage",
        "phage_plasmid",
    ]

    boundaries = ["terminal_inverted_repeat_element", "direct_repeat_element"]

    with open(mobannot) as input_file:
        for line in input_file:
            line = line.rstrip()
            # Annotation lines have exactly 9 columns
            if len(line.split("\t")) == 9:
                data_list = gff_parser(line)
                contig = data_list[0]
                start = data_list[2]
                end = data_list[3]
                if data_list[1] in mges:
                    if data_list[1] == "conjugative_integron":
                        data_list[1] = "conjugative_element"
                    if data_list[1] == "viral_sequence":
                        data_list[1] = "phage"
                    if data_list[1] == "prophage":
                        data_list[1] = "phage"
                    data_list = tuple(data_list)
                    if contig in momofy_dict:
                        momofy_dict[contig].append(data_list)
                    else:
                        momofy_dict[contig] = [data_list]

                    for attribute in line.split("\t")[8].split(";"):
                        att_key, att_val = attribute.split("=")
                        if att_key == "mobile_element_type":
                            min_key = (
                                contig
                                + ":"
                                + data_list[4].split("-")[1].replace(":", "-")
                            )
                            momo_subtypes[min_key] = att_val
                            bm1_key = (contig, start)
                            bound_map_1[bm1_key] = att_val
                            bm2_key = (contig, end)
                            bound_map_2[bm2_key] = att_val

                elif data_list[1] in boundaries:
                    for attribute in line.split("\t")[8].split(";"):
                        att_key, att_val = attribute.split("=")
                        if att_key == "flanking_site":
                            location = att_val.split("_")[1]
                            if location == "1":
                                loc_key = (contig, start)
                                data_list.append(att_val)
                            elif location == "2":
                                loc_key = (contig, end)
                                data_list.append(att_val)
                    irdr_map[loc_key] = data_list

    # Generating the irdr list with the right key: DR_1.1(IME)
    boundaries_counter = 0
    for contig in momofy_dict:
        for element in momofy_dict[contig]:
            boundaries_counter += 1
            start = element[2]
            end = element[3]
            bm1_key = (contig, start)
            bm2_key = (contig, end)
            if all([bm1_key in bound_map_1, bm1_key in irdr_map]):
                new_id = (
                    irdr_map[bm1_key][5]
                    + "."
                    + str(boundaries_counter)
                    + "("
                    + bound_map_1[bm1_key]
                    + ")"
                )
                min_info = [
                    irdr_map[bm1_key][0],
                    irdr_map[bm1_key][1],
                    irdr_map[bm1_key][2],
                    irdr_map[bm1_key][3],
                ]
                irdr_dict[new_id] = min_info
            if all([bm2_key in bound_map_2, bm2_key in irdr_map]):
                new_id = (
                    irdr_map[bm2_key][5]
                    + "."
                    + str(boundaries_counter)
                    + "("
                    + bound_map_2[bm2_key]
                    + ")"
                )
                min_info = [
                    irdr_map[bm2_key][0],
                    irdr_map[bm2_key][1],
                    irdr_map[bm2_key][2],
                    irdr_map[bm2_key][3],
                ]
                irdr_dict[new_id] = min_info

    return (momofy_dict, irdr_dict, momo_subtypes)


def promge_parser(promge, meta):
    mge_type = [
        ("is_tn", "insertion_sequence"),
        ("phage", "phage"),
        ("phage_like", "phage"),
        ("ce", "conjugative_element"),
        ("integron", "integron"),
        ("mi", "mobility_island "),
        ("cellular", "cellular_recombinase"),
    ]
    mge_desc = {}
    with open(meta) as input_meta:
        next(input_meta)
        for line in input_meta:
            (
                is_tn,
                phage,
                phage_like,
                ce,
                integron,
                mi,
                cellular,
                contig,
                start,
                end,
                size,
                n_genes,
                mgeR,
            ) = line.rstrip().split("\t")
            new_id = contig + ":" + start + "-" + end
            description = []
            for index in range(7):
                if int(line.rstrip().split("\t")[index]) > 0:
                    description.append(mge_type[index][1])
            description = "|".join(description)
            mge_desc[new_id] = description

    promge_dict = {}
    mgeR = {}
    with open(promge) as input_gff:
        for line in input_gff:
            line = line.rstrip()
            # Annotation lines have exactly 9 columns
            if len(line.split("\t")) == 9:
                data_list = gff_parser(line)
                contig = data_list.pop(0)
                seq_id = data_list.pop(-1)
                new_id = re.sub(".*contig_", "contig_", seq_id)
                if new_id in mge_desc:
                    desc = mge_desc[new_id]
                    data_list.pop(0)
                    data_list.insert(0, desc)
                else:
                    print("No description for " + seq_id)
                data_list.append(new_id)
                data_list.insert(0, contig)
                data_list = tuple(data_list)
                if contig in promge_dict:
                    promge_dict[contig].append(data_list)
                else:
                    promge_dict[contig] = [data_list]

                for attribute in line.split("\t")[8].split(";"):
                    att_key, att_val = attribute.split("=")
                    if att_key == "mgeR":
                        mgeR[new_id] = att_val

    return (promge_dict, mgeR)


def mapper(momofy_dict, promge_dict):
    momo_unique = []
    pro_unique = []
    momo_used = []
    pro_used = []
    overlapped = []
    permge_metadata = {}

    pro_all = []
    for contig in promge_dict:
        for pro_element in promge_dict[contig]:
            pro_seq_id = pro_element[4]
            pro_all.append(pro_seq_id)
            permge_metadata[pro_seq_id] = pro_element
            # print('proMGE',pro_seq_id, pro_element)

    momo_all = []
    for contig in momofy_dict:
        for momo_element in momofy_dict[contig]:
            momo_coords = momo_element[4].split("-")[1].replace(":", "-")
            momo_seq_id = contig + ":" + momo_coords
            momo_all.append(momo_seq_id)
            permge_metadata[momo_seq_id] = momo_element
            # print('momo',momo_seq_id, momo_element)

    contigs_list = list(set(list(momofy_dict.keys()) + list(promge_dict.keys())))
    for contig in contigs_list:
        if contig in promge_dict:
            for pro_element in promge_dict[contig]:
                pro_seq_type = pro_element[1]
                pro_start = pro_element[2]
                pro_end = pro_element[3]
                pro_seq_id = pro_element[4]
                pro_len = int(pro_end) - int(pro_start)
                pro_range = range(int(pro_start), int(pro_end) + 1)
                flag = 0

                if pro_seq_id not in pro_used:
                    if contig in momofy_dict:
                        for momo_element in momofy_dict[contig]:
                            momo_seq_type = momo_element[1]
                            momo_start = momo_element[2]
                            momo_end = momo_element[3]
                            momo_seq_id = contig + ":" + momo_start + "-" + momo_end
                            momo_len = int(momo_end) - int(momo_start)
                            momo_range = range(int(momo_start), int(momo_end) + 1)
                            intersection = len(list(set(pro_range) & set(momo_range)))
                            if intersection > 0:
                                flag = 1
                                momo_used.append(momo_seq_id)
                                pro_used.append(pro_seq_id)
                                pair = [pro_seq_id, momo_seq_id]
                                cluster_flag = 0
                                index = -1
                                for cluster in overlapped:
                                    index += 1
                                    if any(
                                        [pro_seq_id in cluster, momo_seq_id in cluster]
                                    ):
                                        cluster_flag = 1
                                        overlapped[index].append(pro_seq_id)
                                        overlapped[index].append(momo_seq_id)

                                if cluster_flag == 0:
                                    overlapped.append(pair)

                if flag == 0:
                    pro_unique.append(pro_seq_id)

        momo_unique = list(set(momo_all) - set(momo_used))

    ## Collapsing overlapping clusters
    # saving cluster elements coordinates
    cluster_counter = 0
    final_overlapped = []
    clusters_dict = {}
    clusters_labels = {}
    for cluster in overlapped:
        cluster_counter += 1
        cluster = list(set(cluster))
        # print('clstr_'+str(cluster_counter)+' '+','.join(cluster)+'\n')
        for element in cluster:
            contig = element.split(":")[0]
            start = int(element.split(":")[1].split("-")[0])
            end = int(element.split(":")[1].split("-")[1])
            composite_key = (contig, cluster_counter)
            # print(contig+', clstr_'+str(cluster_counter)+', '+str(start)+', '+str(end))
            if composite_key in clusters_dict:
                clusters_dict[composite_key].append(start)
                clusters_dict[composite_key].append(end)
            else:
                clusters_dict[composite_key] = [start, end]
        clusters_labels[composite_key] = cluster

    # saving clusters boundaries
    clust_bounderies = {}
    for comp_key in clusters_dict:
        min_coord = sorted(clusters_dict[comp_key])[0]
        max_coord = sorted(clusters_dict[comp_key])[-1]
        coords = (min_coord, max_coord)
        clust_bounderies[comp_key] = coords
        # print(comp_key,min_coord,max_coord)

    # Finding overlapping clusters
    clst_used = []
    clst_2_add = {}
    for comp_key_1 in clust_bounderies:
        contig_1 = comp_key_1[0]
        cluster_1 = comp_key_1[1]
        clst_1_start = clust_bounderies[comp_key_1][0]
        clst_1_end = clust_bounderies[comp_key_1][1]
        clst_1_range = (int(clst_1_start), int(clst_1_end) + 1)
        for comp_key_2 in clust_bounderies:
            contig_2 = comp_key_2[0]
            cluster_2 = comp_key_2[1]
            if contig_1 == contig_2:
                if comp_key_1 != comp_key_2:
                    clst_2_start = clust_bounderies[comp_key_2][0]
                    clst_2_end = clust_bounderies[comp_key_2][1]
                    clst_2_range = (int(clst_2_start), int(clst_2_end) + 1)
                    intersection = len(list(set(clst_1_range) & set(clst_2_range)))
                    if intersection > 0:
                        if all(
                            [comp_key_1 not in clst_used, comp_key_2 not in clst_used]
                        ):
                            boundaries_list = [
                                clst_1_start,
                                clst_1_end,
                                clst_2_start,
                                clst_2_end,
                            ]
                            new_min = sorted(boundaries_list)[0]
                            new_max = sorted(boundaries_list)[1]
                            collapsed_cluster = list(
                                set(
                                    clusters_labels[comp_key_1]
                                    + clusters_labels[comp_key_2]
                                )
                            )
                            cluster_counter += 1
                            collapsed_key = (contig_1, cluster_counter)
                            clst_2_add[collapsed_key] = collapsed_cluster
                            # print(comp_key_1, comp_key_2, new_min, new_max, collapsed_cluster)
                            clst_used.append(comp_key_1)
                            clst_used.append(comp_key_2)

    for comp_key in clst_used:
        del clusters_labels[comp_key]

    for comp_key in clst_2_add:
        clusters_labels[comp_key] = clst_2_add[comp_key]

    for cluster in clusters_labels:
        final_overlapped.append(clusters_labels[cluster])

    return (momo_unique, pro_unique, final_overlapped, permge_metadata)


def to_print(metadata_tuple, genome, source, extra_attributes):
    contig = metadata_tuple[0]
    mge_type = metadata_tuple[1].replace("|", "/")
    start = metadata_tuple[2]
    end = metadata_tuple[3]
    element_id = metadata_tuple[4]
    global_id = "ID=" + genome + "|" + contig + ":" + start + "-" + end
    extra_attributes.insert(0, global_id)
    attributes = ";".join(extra_attributes)
    line = (
        "\t".join(
            [
                contig,
                source,
                mge_type,
                start,
                end,
                ".",
                ".",
                ".",
                attributes,
            ]
        )
        + "\n"
    )

    return line


def merger(
    momo_unique,
    pro_unique,
    final_overlapped,
    permge_metadata,
    irdr_dict,
    momo_subtypes,
    mgeR,
    genome_name,
):
    with open(genome_name + "_merged.gff", "w") as to_merged:
        to_merged.write("##gff-version 3\n")
        for mge in pro_unique:
            source = "promge"
            extra_attributes = [
                "merged_label=promge_unique",
                "merged_coords=NA",
                "merged_types=NA",
                "subtype=NA",
                "mge_recombinase=" + mgeR[mge],
            ]
            line = to_print(permge_metadata[mge], genome_name, source, extra_attributes)
            to_merged.write(line)

        for mge in momo_unique:
            source = "MAP"
            extra_attributes = [
                "merged_label=MAP_unique",
                "merged_coords=NA",
                "merged_types=NA",
                "subtype=" + momo_subtypes[mge],
                "mge_recombinase=NA",
            ]
            line = to_print(permge_metadata[mge], genome_name, source, extra_attributes)
            to_merged.write(line)

        for boundary in irdr_dict:
            source = "MAP"
            extra_attributes = [
                "merged_label=MAP_unique",
                "merged_coords=NA",
                "merged_types=NA",
                "subtype=NA",
                "mge_recombinase=NA",
            ]

            contig = irdr_dict[boundary][0]
            mge_type = irdr_dict[boundary][1]
            start = irdr_dict[boundary][2]
            end = irdr_dict[boundary][3]
            global_id = "ID=" + boundary
            extra_attributes.insert(0, global_id)
            attributes = ";".join(extra_attributes)
            line = (
                "\t".join(
                    [
                        contig,
                        source,
                        mge_type,
                        start,
                        end,
                        ".",
                        ".",
                        ".",
                        attributes,
                    ]
                )
                + "\n"
            )
            to_merged.write(line)

        for cluster in final_overlapped:
            cluster = list(set(cluster))
            starts = []
            ends = []
            coord_list = []
            nested_types = []
            mgeR_list = []
            momosub_list = []
            promge_positions = []
            momofy_positions = []

            for element in cluster:
                contig = element.split(":")[0]
                element_type = permge_metadata[element][1].replace(" ", "")
                if "|" in element_type:
                    for composite_type in element_type.split("|"):
                        nested_types.append(composite_type)
                else:
                    nested_types.append(element_type)

                mge_start = int(element.split(":")[1].split("-")[0])
                starts.append(mge_start)
                mge_end = int(element.split(":")[1].split("-")[1])
                ends.append(mge_end)

                element_range = list(range(int(mge_start), int(mge_end) + 1))
                if element in mgeR:
                    mgeR_list.append(mgeR[element])
                    promge_positions = promge_positions + element_range

                if element in momo_subtypes:
                    momosub_list.append(momo_subtypes[element])
                    momofy_positions = momofy_positions + element_range

                mge_data = element_type + ":" + element.split(":")[1]
                coord_list.append(mge_data)

            source = "promge/MAP"
            start = sorted(starts)[0]
            end = sorted(ends)[-1]

            # Finding coverage per method
            promge_positions = len(list(set(promge_positions)))
            momofy_positions = len(list(set(momofy_positions)))
            merged_len = int(end) - int(start)
            promge_cov = float(promge_positions) / float(merged_len)
            momo_cov = float(momofy_positions) / float(merged_len)

            if all([promge_cov >= 0.75, momo_cov >= 0.75]):
                mge_type = "complete_overlap"
            else:
                mge_type = "partial_overlap"

            global_id = (
                "ID=" + genome_name + "|" + contig + ":" + str(start) + "-" + str(end)
            )

            merge_info = ",".join(coord_list)
            nested_info = ",".join(list(set(nested_types)))
            mgeR_info = ",".join(mgeR_list)
            momo_sub_info = ",".join(momosub_list)

            extra_attributes = [
                "merged_label=" + mge_type,
                "merged_coords=" + merge_info,
                "merged_types=" + nested_info,
                "subtype=" + momo_sub_info,
                "mge_recombinase=" + mgeR_info,
            ]

            extra_attributes.insert(0, global_id)
            attributes = ";".join(extra_attributes)
            line = (
                "\t".join(
                    [
                        contig,
                        source,
                        "nested",
                        str(start),
                        str(end),
                        ".",
                        ".",
                        ".",
                        attributes,
                    ]
                )
                + "\n"
            )

            to_merged.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="This script merge the coordinates of predicted MGEs using the Mobilome Annotation Pipeline v2.0 and proMGE v2023. Please provide the gff files containing the coordinates of each prediction and the metadata file for proMGE."
    )
    parser.add_argument(
        "--mobannot",
        type=str,
        help="Mobilome Annotation Pipeline v2.0 predictions (no genes) (gff)",
        required=True,
    )
    parser.add_argument(
        "--proMGE",
        type=str,
        help="ProMGE annotation file (gff)",
        required=True,
    )
    parser.add_argument(
        "--meta",
        type=str,
        help="ProMGE metadata file (tsv)",
        required=True,
    )
    parser.add_argument(
        "--genome_name",
        type=str,
        help="Genome name will be used to build the mge id and as prefix for the output file",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    (momofy_dict, irdr_dict, momo_subtypes) = momo_parser(args.mobannot)
    (promge_dict, mgeR) = promge_parser(args.proMGE, args.meta)

    (momo_unique, pro_unique, final_overlapped, permge_metadata) = mapper(
        momofy_dict, promge_dict
    )

    merger(
        momo_unique,
        pro_unique,
        final_overlapped,
        permge_metadata,
        irdr_dict,
        momo_subtypes,
        mgeR,
        args.genome_name,
    )


if __name__ == "__main__":
    main()
