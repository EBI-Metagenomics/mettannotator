#!/usr/bin/env python3

import argparse
import logging
import os
import sys

logging.basicConfig(level=logging.INFO)


def main(infile, outdir):
    taxid = assign_taxid(infile)
    check_dir(outdir)
    outfile = "proteins.fasta"
    outpath = os.path.join(outdir, outfile)
    with open(outpath, "w") as file_out, open(infile, "r") as file_in:
        for line in file_in:
            if line.startswith(">"):
                formatted_line = reformat_line(line, taxid)
                file_out.write(formatted_line)
            else:
                file_out.write(line)


def check_dir(directory_path):
    if not os.path.exists(directory_path):
        try:
            os.makedirs(directory_path)
        except OSError as e:
            logging.error(f"Error: Failed to create directory '{directory_path}'. {e}")


def reformat_line(line, taxid):
    line = line.lstrip('>').strip()
    id, description = line.split(maxsplit=1)
    if taxid == "820":
        sp_name = "Bacteroides uniformis"
    elif taxid == "821":
        sp_name = "Phocaeicola vulgatus"
    elif taxid == "46503":
        sp_name = "Parabacteroides merdae"
    else:
        raise ValueError("Unknown species")
    formatted_line = ">tr|{id}|{description} OS={sp_name} OX={taxid}\n".format(id=id, description=description,
                                                                               sp_name=sp_name, taxid=taxid)
    return formatted_line


def assign_taxid(infile):
    try:
        with open(infile, 'r') as file:
            # Read the first line
            first_line = file.readline().strip()
            species_code = first_line[1:3]

            # Assign taxid based on species code
            if species_code == "BU":
                taxid = "820"
            elif species_code == "PV":
                taxid = "821"
            elif species_code == "PM":
                taxid = "46503"
            else:
                raise ValueError("Unknown species")
            return taxid
    except Exception as e:
        logging.error(f"Error: {e}")
        exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script reformats the fasta faa file to prepare it for UniRule."
        )
    )
    parser.add_argument(
        "-i",
        dest="infile",
        required=True,
        help="Input protein fasta file.",
    )
    parser.add_argument(
        "-o",
        dest="outdir",
        required=True,
        help=(
            "Path to the folder where the output will be saved to."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.infile,
        args.outdir,
    )
