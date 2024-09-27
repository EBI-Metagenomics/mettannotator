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
import os
import sys

import requests
from retry import retry

logging.basicConfig(level=logging.INFO)


def main(taxid, include_kingdom, outfile):
    kingdom = "Bacteria"
    if not taxid.isdigit():
        sys.exit(f"Invalid format of taxid {taxid}")
    try:
        url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}"
        r = run_request(url)
        res = r.json()
        lineage = res.get("lineage", "")
        rank = res.get("rank", "")
        if rank == "superkingdom":
            lineage = res.get("scientificName", "")
        if lineage.startswith("Archaea"):
            kingdom = "Archaea"
        elif lineage.startswith("Vir"):
            kingdom = "Viruses"
        elif lineage.startswith("Bac"):
            pass
        else:
            logging.error(
                f"Unknown lineage {lineage}. Reporting default kingdom instead: Bacteria."
            )
    except:  # noqa E722
        logging.error(
            f"Unable to identify lineage for taxid {taxid}. "
            "Reporting default kingdom instead: Bacteria."
        )
    logging.info(f"Reporting kingdom {kingdom} for taxid {taxid}")
    if outfile:
        if include_kingdom:
            outfile_path, outfile_name = os.path.split(outfile)
            outfile_new_name = f"{kingdom}_{outfile_name}"
            outfile = os.path.join(outfile_path, outfile_new_name)
        with open(outfile, "w") as file_out:
            file_out.write(kingdom)
    else:
        print(kingdom.strip())


@retry(tries=5, delay=10, backoff=1.5)
def run_request(full_url):
    r = requests.get(url=full_url)
    r.raise_for_status()
    return r


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes taxid as input and returns the kingdom. If kingdom cannot be inferred "
            "from the taxid, the script returns 'Bacteria' by default. Only works for prokaryotic and "
            "viral genomes."
        )
    )
    parser.add_argument(
        "-t",
        dest="taxid",
        required=True,
        help="Taxid to identify the kingdom for.",
    )
    parser.add_argument(
        "--include-kingdom",
        action="store_true",
        help="Add the kingdom name at the beginning of the file name.",
    )
    parser.add_argument(
        "-o",
        dest="outfile",
        required=False,
        help=(
            "Path to the output file (with file name) where the result will be saved. If the --include-kingdom option "
            "was specified, the kingdom will be added to the beginning of the chosen output file name. "
            "If the output file is not specified, the script will print the output as STDOUT."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.taxid,
        args.include_kingdom,
        args.outfile,
    )
