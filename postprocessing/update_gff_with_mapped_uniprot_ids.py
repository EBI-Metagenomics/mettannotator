#!/usr/bin/env python3

import sys
import os
import subprocess


def run_blast(pipeline_fasta, uniprot_fasta, blast_output):
    # Create blast database
    subprocess.run(
        ["makeblastdb", "-in", pipeline_fasta, "-parse_seqids", "-dbtype", "prot"]
    )

    # Run blastp
    subprocess.run(
        [
            "blastp",
            "-db",
            pipeline_fasta,
            "-query",
            uniprot_fasta,
            "-out",
            blast_output,
            "-outfmt",
            "6",
            "-evalue",
            "1e-10",
        ]
    )


def extract_best_hits(blast_output):
    best_hits = {}
    with open(blast_output, "r") as file:
        for line in file:
            fields = line.strip().split("\t")
            query_id = fields[0]
            result_id = fields[1]
            percent_id = fields[2]
            overlap = fields[3]
            if query_id not in best_hits or float(percent_id) > float(
                best_hits[query_id][2]
            ):
                best_hits[query_id] = (result_id, percent_id, overlap)
    return best_hits


def create_mapping_file(best_hits, mapping_file):
    with open(mapping_file, "w") as file:
        file.write("Query\tResult\tPercent_ID\tOverlap\n")
        for query_id, values in best_hits.items():
            result_id, percent_id, overlap = values
            file.write(f"{query_id}\t{result_id}\t{percent_id}\t{overlap}\n")


def update_gff_with_mapping(gff_file, mapping_file, output_file):
    mapping_dict = {}
    with open(mapping_file, "r") as file:
        next(file)  # Skip header
        for line in file:
            query, result, _, _ = line.strip().split("\t")
            mapping_dict[result] = query

    with open(gff_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fields = line.strip().split("\t")
                if (
                    len(fields) >= 9
                    and fields[2] == "CDS"
                    and fields[8].startswith("ID=")
                ):
                    # Should do something better than assume ID is always first
                    identifier = fields[8].split(";")[0].split("=")[1]
                    if identifier.startswith("CDS:"):
                        identifier = identifier[len("CDS:") :]
                    if identifier in mapping_dict:
                        fields[8] += f";Dbxref=UniProt:{mapping_dict[identifier]}"
                outfile.write("\t".join(fields) + "\n")


if __name__ == "__main__":

    if len(sys.argv) != 5:
        print(
            """
                Usage: python update_gff_with_mapped_uniprot_ids.py <species> <gff_file> <uniprot.faa> <pipeline.faa>

		Need to download proteomes from UniProt (Buniformis: UP000004110, Pvulgatus: UP000002861) & strip away any unwanted bits in the headers, and have the new GFF and protein fasta files.
		Blastp runs on the two protein fasta files, with evalue -10. The best hits are extracted from the results to produce a mapping file between the new genes and the old uniprot identifiers.
		This mapping is used to insert Dbxref=UniProt: key-value attributes to the ninth column of CDS lines in the final output.
                The species name is just used to name files.
        """
        )
        sys.exit(1)

    species = sys.argv[1]
    gff_file = sys.argv[2]
    uniprot_fasta = sys.argv[3]
    pipeline_fasta = sys.argv[4]
    blast_output = species + ".blastp.evalue.out"
    mapping_file = species + ".mapping.txt"
    output_gff_file = species + ".withuniprotids.output.gff"

    # Run blast if blast results do not already exist
    if not os.path.exists(blast_output):
        run_blast(pipeline_fasta, uniprot_fasta, blast_output)

    # Extract best hits
    best_hits = extract_best_hits(blast_output)

    # Create mapping file
    create_mapping_file(best_hits, mapping_file)

    # Update GFF with mapping
    update_gff_with_mapping(gff_file, mapping_file, output_gff_file)
