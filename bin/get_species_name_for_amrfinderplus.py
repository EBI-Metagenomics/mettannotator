#!/usr/bin/env python3

# Copyright 2023-2025 EMBL - European Bioinformatics Institute
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
import os

from requests import Response
import requests
from retry import retry
from dataclasses import dataclass


def main(taxid: str, database_folder: str):
    # Load taxonomies that AMRFinderPlus accepts
    amrfinderplus_species_list, amrfinderplus_genus_list = load_organisms(database_folder)

    # Check the taxonomy of the provided taxid
    taxonomy = fetch_taxonomy(taxid)

    organism = deduce_organism(taxonomy, amrfinderplus_species_list, amrfinderplus_genus_list)
    print(organism or "")  # print empty string if no match


@dataclass
class Taxonomy:
    rank: str
    scientific_name: str
    lineage: str


def fetch_taxonomy(taxid: str) -> Taxonomy:
    """Fetch taxonomy info from ENA API for a given taxid."""
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxid}"
    r = run_request(url)
    res = r.json()
    rank = res.get("rank", "")
    scientific_name = res.get("scientificName", "")
    lineage = res.get("lineage", "")

    if rank in {"no rank", "strain"} and lineage:
        # take the next lowest taxon in lineage to see if we can get to a species or a genus
        next_lowest_taxon = [p.strip() for p in lineage.split(";") if p.strip()][-1]
        name_url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/scientific-name/{next_lowest_taxon}"
        r2 = run_request(name_url)
        res2 = r2.json()[0]
        rank = res2.get("rank", "")
        scientific_name = res2.get("scientificName", "")
        lineage = res2.get("lineage", "")
    return Taxonomy(rank=rank, scientific_name=scientific_name, lineage=lineage)


def deduce_organism(taxonomy: Taxonomy, species_list: list, genus_list: list):
    """Deduce AMRFinderPlus-compatible organism name from taxonomy.

    Returns a string if matched, otherwise None.
    """
    name = taxonomy.scientific_name.strip()

    # Report Shigella → Escherichia per AMRFinderPlus instructions
    if name.lower().startswith(("shigella", "escherichia")) and taxonomy.rank in ["genus", "species"]:
        return "Escherichia"

    if taxonomy.rank == "genus":
        return name if name in genus_list else None

    if taxonomy.rank == "species":
        merged_name = name.replace(" ", "_")
        return merged_name if merged_name in species_list else None

    return None


@retry(tries=5, delay=10, backoff=1.5)
def run_request(full_url: str) -> Response:
    r = requests.get(url=full_url)
    r.raise_for_status()
    return r


def load_organisms(database_folder) -> [list, list]:
    """Deduce AMRFinderPlus-compatible organism name from taxonomy."""
    amrfinderplus_species_list = list()
    amrfinderplus_genus_list = list()
    tax_file = os.path.join(database_folder, "taxgroup.tsv")
    if not os.path.isfile(tax_file):
        raise FileNotFoundError(f"Cannot find the expected tax file in AMRFinder db: {tax_file}")
    with open(tax_file, "r") as file_in:
        for line in file_in:
            if line.startswith("#"):
                continue
            taxon_name = line.split("\t")[0]
            # check if this is a species or a genus
            if "_" in taxon_name:
                amrfinderplus_species_list.append(taxon_name)
            else:
                amrfinderplus_genus_list.append(taxon_name)
    return amrfinderplus_species_list, amrfinderplus_genus_list


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "The script takes taxid as input and returns the organism name that can be used for AMRFinderPlus. If "
            "no organism name is available (taxonomy is not species-level or species is not in the AMRFinderPlus "
            "list, the script returns an empty string."
        )
    )
    parser.add_argument(
        "-t",
        dest="taxid",
        required=True,
        help="Taxid to identify the organism for.",
    )
    parser.add_argument(
        "-d",
        dest="database_folder",
        required=False,
        help=(
            "Path to the AMRFinderPlus database folder."
        ),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(
        args.taxid,
        args.database_folder,
    )
