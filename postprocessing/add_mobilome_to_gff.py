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


def main(mobilome_file, infile, outfile):
    mobilome_dict = load_mobilome(mobilome_file)
    printed_list = (
        list()
    )  # record which mobilome lines have already been printed to the output file
    fasta_flag = False
    with open(infile) as file_in, open(outfile, "w") as file_out:
        previous_contig = ""
        for line in file_in:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    fasta_flag = True
                file_out.write(line)
            elif fasta_flag:
                # We are printing the FASTA sequence at the end of the file now
                file_out.write(line)
            else:
                contig, tool, feature, start, end, blank1, blank2, blank3, col9 = (
                    line.strip().split("\t")
                )
                if previous_contig and contig != previous_contig:
                    # We switched to a new contig, check if there are any MGE lines left on the previous contig that
                    # were not printed.
                    last_contig_lines_to_print = look_for_lines_to_print(
                        mobilome_dict, previous_contig, 9999999999999, printed_list
                    )
                    for line_to_print in last_contig_lines_to_print:
                        file_out.write(line_to_print)
                        printed_list.append(line_to_print)
                previous_contig = contig

                # Now process the current line, check if there is an overlap with any MGEs
                start, end = int(start), int(end)
                overlap, lines_to_print = check_overlap(
                    mobilome_dict, contig, start, end
                )
                # Before printing any results, check if there are earlier mobilome lines on this contig
                # that haven't been printed yet
                extra_lines_to_print = look_for_lines_to_print(
                    mobilome_dict, contig, start, printed_list
                )
                for extra_line in extra_lines_to_print:
                    if extra_line not in lines_to_print:
                        file_out.write(extra_line)
                        printed_list.append(extra_line)
                if overlap:
                    add_to_cds = ""
                    mge_ids = list()
                    mge_types = list()
                    # Add mobilome record info to CDS
                    for line_to_print in lines_to_print:
                        if line_to_print not in printed_list:
                            # print the mobilome line first
                            file_out.write(line_to_print)
                            printed_list.append(line_to_print)
                        if feature == "CDS":
                            # extract information to add to the CDS
                            mge_col9 = line_to_print.strip().split("\t")[8]
                            attributes_dict = dict(
                                re.split(r"(?<!\\)=", item)
                                for item in re.split(r"(?<!\\);", mge_col9)
                            )
                            mge_ids.append(attributes_dict["ID"])
                            if "merged_types" in attributes_dict:
                                mge_types.append(attributes_dict["merged_types"])
                            else:
                                mge_types.append(line_to_print.strip().split("\t")[2])
                    if feature == "CDS":
                        if mge_ids:
                            add_to_cds += "mge_id={}".format(",".join(mge_ids))
                        if mge_types:
                            add_to_cds += ";mge_types={}".format(",".join(mge_types))
                        if add_to_cds and not col9.endswith(";"):
                            col9 += ";"
                        col9 = col9 + add_to_cds
                    line = (
                        "\t".join(
                            [
                                contig,
                                tool,
                                feature,
                                str(start),
                                str(end),
                                blank1,
                                blank2,
                                blank3,
                                col9,
                            ]
                        )
                        + "\n"
                    )
                    file_out.write(line)
                else:
                    file_out.write(line)
    sanity_check(mobilome_dict, printed_list)


def sanity_check(mobilome_dict, printed_list):
    """Check that the number of records added to the GFF matches the number of records in the mobilome file"""
    printed_list_length = len(list(set(printed_list)))
    mobilome_count = sum(len(inner_dict) for inner_dict in mobilome_dict.values())
    assert (
        printed_list_length == mobilome_count
    ), f"The number of mobilome entries added to the GFF does not match the expected count: added {printed_list_length}, expected{mobilome_count}"


def look_for_lines_to_print(mobilome_dict, contig, start, printed_list):
    lines_to_print = list()
    if contig not in mobilome_dict:
        return lines_to_print
    for interval in mobilome_dict[contig]:
        if interval[0] < start and mobilome_dict[contig][interval] not in printed_list:
            lines_to_print.append(mobilome_dict[contig][interval])
    return lines_to_print


def calculate_overlap_fraction(interval1, interval2):
    """Calculate what percentage of interval 1 is contained inside interval 2"""
    start1, end1 = interval1
    start2, end2 = interval2

    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)

    if overlap_start < overlap_end:
        overlap_length = overlap_end - overlap_start + 1
        interval1_length = end1 - start1 + 1
        overlap_fraction = overlap_length / interval1_length
        return overlap_fraction
    else:
        return 0.0


def check_overlap(dictionary, sequence, start, end):
    """Check if at least 75% of the CDS overlaps with any entry in the MGE dictionary."""
    result = list()
    if sequence in dictionary:
        for region, text in dictionary[sequence].items():
            if calculate_overlap_fraction((start, end), region) >= 0.75:
                result.append(text)
    if len(result) > 0:
        return True, result
    else:
        return False, ""


def load_mobilome(infile):
    mobilome_dict = dict()
    with open(infile) as file_in:
        for line in file_in:
            if line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                else:
                    continue
            contig, tool, feature, start, end, blank1, blank2, blank3, col9 = (
                line.strip().split("\t")
            )
            if feature == "nested":
                # replace with a more meaningful name
                feature = "nested_mobile_element"
            else:
                # remove merged attributes from col9 of records that are not merged
                col9 = (
                    col9.replace("merged_coords=NA;", "")
                    .replace("merged_types=NA;", "")
                    .replace("subtype=NA;", "")
                )
                if feature in [
                    "terminal_inverted_repeat_element",
                    "direct_repeat_element",
                ]:
                    col9 = col9.replace(";mge_recombinase=NA", "")
            line = (
                "\t".join(
                    [contig, tool, feature, start, end, blank1, blank2, blank3, col9]
                )
                + "\n"
            )
            mobilome_dict.setdefault(contig, dict()).update(
                {tuple([int(start), int(end)]): line}
            )
    return mobilome_dict


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes the mobilome GFF and the genome annotation GFF produced by mettannotator "
            "as inputs and adds mobilome to the main GFF. The current version of the script is project specific "
            "and is not meant to be a part of the general annotation workflow."
        )
    )
    parser.add_argument(
        "-m",
        dest="mobilome",
        required=True,
        help="GFF containing the mobilome records.",
    )
    parser.add_argument(
        "-i",
        dest="infile",
        required=True,
        help="The GFF annotation file produced by mettannotator.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help="The name of the file the result will be saved to.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.mobilome,
        args.infile,
        args.outfile,
    )
