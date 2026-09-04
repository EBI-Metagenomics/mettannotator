#!/usr/bin/env python3
# Normalise an Ensembl Bacteria GFF3 so that annotate_gff.py can process it.
#
# Ensembl GFF3 problems this script resolves:
#   - CDS ID is "CDS:<protein_id>" — strip to bare protein_id
#   - CDS features lack locus_tag and product — pulled from parent gene feature
#   - ### separator lines accumulate as false header lines in annotate_gff.py
#   - Ensembl-specific pragma headers (#!) are not informative for mettannotator
#
# Output: standard Prokka-style GFF3 with ID, locus_tag, product on every CDS.

import argparse
import sys
import urllib.parse

SKIP_FEATURES = {
    "mRNA",
    "exon",
    "pseudogene",
    "pseudogenic_transcript",
    "ncRNA_gene",
    "ncRNA",
    "pre_miRNA",
    "miRNA",
    "snoRNA",
    "snRNA",
    "lnc_RNA",
    "antisense_RNA",
    "ribozyme",
    "rRNA",
    "tRNA",
    "biological_region",
    "chromosome",
    "region",
    "tmRNA",
}


def parse_col9(col9):
    """Return an ordered list of (key, value) pairs from a GFF3 col9 string."""
    pairs = []
    for item in col9.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            pairs.append((k, v))
        else:
            pairs.append((item, ""))
    return pairs


def col9_to_dict(col9):
    return {k: v for k, v in parse_col9(col9)}


def decode_description(text):
    """URL-decode an Ensembl description field and sanitise for GFF3 col9."""
    text = urllib.parse.unquote(text)
    # Semicolons are col9 attribute separators — replace to avoid parse errors
    text = text.replace(";", ",")
    return text.strip()


def normalise(in_gff, out_gff):
    gene_info = {}  # gene_id  -> {"name": str|None, "description": str|None}
    mrna_to_gene = {}  # transcript_id -> gene_id

    # --- Pass 1: collect gene and mRNA metadata ---
    with open(in_gff) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line:
                continue
            cols = line.split("\t")
            if len(cols) != 9:
                continue
            feature = cols[2]
            attrs = col9_to_dict(cols[8])

            if feature == "gene":
                # gene_id attribute is the stable identifier (e.g. ABAYE0001)
                gid = attrs.get("gene_id") or attrs.get("ID", "").removeprefix("gene:")
                desc = decode_description(attrs.get("description", ""))
                name = attrs.get("Name")
                gene_info[gid] = {"name": name, "description": desc or None}

            elif feature == "mRNA":
                tid = attrs.get("transcript_id") or attrs.get("ID", "").removeprefix(
                    "transcript:"
                )
                parent_raw = attrs.get("Parent", "")
                gid = parent_raw.removeprefix("gene:")
                if tid and gid:
                    mrna_to_gene[tid] = gid

    # --- Pass 2: write normalised GFF3 ---
    with open(in_gff) as fh, open(out_gff, "w") as out:
        for line in fh:
            line = line.rstrip("\n")

            if not line:
                continue

            # GFF3 feature separator — discard
            if line == "###":
                continue

            if line.startswith("#"):
                # Keep standard GFF3 directives; drop Ensembl-specific pragmas (#!...)
                if line.startswith("##gff-version") or line.startswith(
                    "##sequence-region"
                ):
                    out.write(line + "\n")
                continue

            cols = line.split("\t")
            if len(cols) != 9:
                continue

            feature = cols[2]

            if feature in SKIP_FEATURES:
                continue

            if feature == "gene":
                attrs = col9_to_dict(cols[8])
                gid = attrs.get("gene_id") or attrs.get("ID", "").removeprefix("gene:")
                name = attrs.get("Name")
                # Match Prokka compliant format: gene ID gets _gene suffix so CDS can
                # reference it via Parent=<gid>_gene
                new_attrs = [("ID", f"{gid}_gene"), ("locus_tag", gid)]
                if name:
                    new_attrs.append(("Name", name))
                cols[8] = ";".join(f"{k}={v}" for k, v in new_attrs)
                out.write("\t".join(cols) + "\n")
                continue

            if feature not in {"CDS", "gene"}:
                print(
                    f"WARNING: Skipping unexpected feature type '{feature}' "
                    f"at {cols[0]}:{cols[3]}-{cols[4]}",
                    file=sys.stderr,
                )
                continue

            attrs = col9_to_dict(cols[8])

            # Stable protein accession — use protein_id if present, else strip CDS: prefix
            protein_id = attrs.get("protein_id")
            if not protein_id:
                raw_id = attrs.get("ID", "")
                protein_id = raw_id.removeprefix("CDS:")

            if not protein_id:
                print(
                    f"WARNING: CDS at {cols[0]}:{cols[3]}-{cols[4]} has no protein_id — skipped",
                    file=sys.stderr,
                )
                continue

            # Resolve gene info via transcript → gene chain
            tid = attrs.get("Parent", "").removeprefix("transcript:")
            gid = mrna_to_gene.get(tid, "")
            gene_meta = gene_info.get(gid, {})

            product = gene_meta.get("description") or "hypothetical protein"
            gene_name = gene_meta.get("name")

            # Build clean Prokka-compliant attributes.
            # No Parent link: BCBio would nest CDS under gene as sub-features, causing
            # gbk_generator to omit CDS from the GBK (it only iterates top-level features).
            # Gene + CDS are linked implicitly by shared locus_tag, matching Prokka's GBK format.
            new_attrs = [
                ("ID", protein_id),
                ("locus_tag", gid if gid else protein_id),
                ("product", product),
            ]
            if gene_name:
                new_attrs.append(("Name", gene_name))

            cols[8] = ";".join(f"{k}={v}" for k, v in new_attrs)
            out.write("\t".join(cols) + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Normalise an Ensembl Bacteria GFF3 for mettannotator."
    )
    parser.add_argument("-i", "--input", required=True, help="Input Ensembl GFF3")
    parser.add_argument("-o", "--output", required=True, help="Output normalised GFF3")
    args = parser.parse_args()
    normalise(args.input, args.output)


if __name__ == "__main__":
    main()
