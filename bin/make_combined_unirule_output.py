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
import logging
import os
import sys

logging.basicConfig(level=logging.INFO)


def main(indir, outfile):
    inputs = {
        "ARBA": "predictions_arba.out",
        "UniRule": "predictions_unirule.out",
        "UniRule-PIRSR": "predictions_unirule-pirsr.out"
    }
    header = "Source\tEvidence\tProteinId\tAnnotationType\tValue\tStart\tEnd\n"
    result = dict()
    for dbname, filename in inputs.items():
        result = load_file(dbname, indir, filename, result)

    with open(outfile, "w") as file_out:
        file_out.write(header)
        for prot, lines in result.items():
            for line in lines:
                file_out.write(line)
                

def load_file(dbname, indir, filename, result):
    if not os.path.exists(os.path.join(indir, filename)):
        logging.error("File {} does not exist in folder {}. Aborting.".format(filename, indir))
        sys.exit()
    with open(os.path.join(indir, filename), "r") as file_in:
        for line in file_in:
            if not line.startswith("Evidence"):
                prot = line.split("\t")[1]
                new_line = "{}\t{}".format(dbname, line)
                result.setdefault(prot, list()).append(new_line)
    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script combines UniFIRE output files into one."
        )
    )
    parser.add_argument(
        "-i",
        dest="infolder",
        required=True,
        help="Folder with uniFIRE outputs.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=True,
        help=(
            "Path to the output file."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.infolder,
        args.outfile,
    )
