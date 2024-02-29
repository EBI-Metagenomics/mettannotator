#!/usr/bin/env python3

import sys

def read_coordinates_file(coord_file):
    coordinates_dict = {}
    with open(coord_file, 'r') as file:
        for line in file:
            gene_name, start, end, essentiality, _ = line.strip().split('\t')
            coordinates_dict[(start, end)] = essentiality
    return coordinates_dict

def update_gff_with_status(gff_file, coordinates_dict):
    updated_lines = []
    matched_coordinates = set()
    with open(gff_file, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            if len(columns) >= 9 and columns[2] == "CDS":
                start, end = columns[3], columns[4]
                essentiality = coordinates_dict.get((start, end))
                if essentiality:
                    columns[8] += f";transit_combined_hmm_gumbel_essentiality={essentiality}"
                    matched_coordinates.add((start,end))
            updated_lines.append('\t'.join(columns) + '\n')
        else:
                updated_lines.append(line)
    return updated_lines, matched_coordinates

def write_updated_gff(output_file, updated_lines):
    with open(output_file, 'w') as file:
        for line in updated_lines:
            file.write(line)

def write_not_matched_file(not_matched_file, coordinates_dict, matched_coordinates):
    not_matched_lines = [f"{start}\t{end}\t{essentiality}\n" for (start, end), (essentiality) in coordinates_dict.items() if (start, end) not in matched_coordinates]
    with open(not_matched_file, 'w') as file:
        file.writelines(not_matched_lines)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("""
                Usage: python script.py <gff_file> <essentiality_file> <output_gff_file> <not_matched_file>

                Match coordinates from essentiality file, produce modified GFF with new key-value pair and a file containing unmatched essentiality lines.
                The key inserted into the 9th column of CDS lines only is transit_combined_hmm_gumbel_essentiality (all lower case to fit with GFF specs)
                Note that the input is a summary of essentiality data (i.e. medium ignored for now)
        """)
        sys.exit(1)

    gff_file = sys.argv[1]
    coordinates_file = sys.argv[2]
    output_file = sys.argv[3]
    not_matched_file = sys.argv[4]

    coordinates_dict = read_coordinates_file(coordinates_file)
    updated_lines, matched_coordinates = update_gff_with_status(gff_file, coordinates_dict)
    write_updated_gff(output_file, updated_lines)
    write_not_matched_file(not_matched_file, coordinates_dict, matched_coordinates)
