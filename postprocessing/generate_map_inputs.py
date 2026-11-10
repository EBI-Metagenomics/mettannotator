#!/usr/bin/env python3
"""
Generate MAP (Microbial Annotation Pipeline) input files from mettanotator results.

Produces:
  1. MAP samplesheet         (sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv)
  2. MAP annotation manifest  (sample,amrfinder_tsv,antismash_gff,gecco_gff,sanntis_gff)
  3. A stdout summary (number of samples, how many have virify results)

Inputs:
  1. mettanotator input samplesheet (prefix,assembly,taxid)
  2. path to virify results         (optional)
  3. path to mettanotator results

Sample names come from column 1 (prefix) of the mettanotator samplesheet and must
match the per-sample result directories under the mettanotator results path. If none
match, the provided paths are assumed wrong and the script exits with an error.
"""

import argparse
import csv
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Per-sample relative paths inside the mettanotator results directory.
# Each is relative to  <mettanotator_results>/<sample>/
# {s} is substituted with the sample name.
# ---------------------------------------------------------------------------
METTANOTATOR_PATHS = {
    "proteins_gff":     "functional_annotation/merged_gff/{s}_annotations.gff",
    "proteins_faa":     "functional_annotation/prokka/{s}.faa",
    "interproscan_tsv": "functional_annotation/interproscan/{s}.IPS.tsv",
    "amrfinder_tsv":    "antimicrobial_resistance/amrfinder_plus/{s}_amrfinderplus.tsv",
    "antismash_gff":    "biosynthetic_gene_clusters/antismash/{s}_antismash.gff",
    "gecco_gff":        "biosynthetic_gene_clusters/gecco/{s}_gecco_clusters.gff",
    "sanntis_gff":      "biosynthetic_gene_clusters/sanntis/{s}_sanntis.gff",
}

# virify gff, relative to <virify_results>/<sample>/
VIRIFY_REL = "08-final/gff/{s}_virify.gff.gz"

MAP_SAMPLESHEET_COLS = [
    "sample", "assembly", "proteins_gff", "proteins_faa",
    "virify_gff", "interproscan_tsv",
]
MANIFEST_COLS = [
    "sample", "amrfinder_tsv", "antismash_gff", "gecco_gff", "sanntis_gff",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate MAP samplesheet and annotation manifest from mettanotator results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input-samplesheet", required=True, type=Path,
        help="mettanotator input samplesheet (CSV with columns: prefix,assembly,taxid)",
    )
    p.add_argument(
        "-m", "--mettanotator-results", required=True, type=Path,
        help="path to the mettanotator results directory",
    )
    p.add_argument(
        "-v", "--virify-results", type=Path, default=None,
        help="path to the virify results directory (optional)",
    )
    p.add_argument(
        "-s", "--out-samplesheet", type=Path, default=Path("map_samplesheet.csv"),
        help="output path for the MAP samplesheet",
    )
    p.add_argument(
        "-a", "--out-manifest", type=Path, default=Path("map_annotation_manifest.csv"),
        help="output path for the MAP annotation manifest",
    )
    p.add_argument(
        "--warn-missing", action="store_true",
        help="print a warning to stderr for every expected mettanotator file that is missing",
    )
    return p.parse_args()


def read_input_samplesheet(path):
    """Return a list of (sample, assembly) tuples from the mettanotator samplesheet."""
    if not path.is_file():
        sys.exit(f"ERROR: input samplesheet not found: {path}")

    rows = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            sys.exit(f"ERROR: input samplesheet is empty: {path}")
        # Be tolerant of header naming: use the first two columns by position
        # but prefer named columns if present.
        fields = [f.strip() for f in reader.fieldnames]
        prefix_col = "prefix" if "prefix" in fields else fields[0]
        assembly_col = "assembly" if "assembly" in fields else fields[1]
        for i, row in enumerate(reader, start=2):
            sample = (row.get(prefix_col) or "").strip()
            assembly = (row.get(assembly_col) or "").strip()
            if not sample:
                print(f"WARNING: skipping row {i}: empty sample name", file=sys.stderr)
                continue
            rows.append((sample, assembly))
    if not rows:
        sys.exit(f"ERROR: no sample rows found in input samplesheet: {path}")
    return rows


def resolve(base, rel, sample, warn_missing, label):
    """Return absolute path string if the file exists, else empty string."""
    candidate = base / sample / rel.format(s=sample)
    if candidate.is_file():
        return str(candidate.resolve())
    if warn_missing:
        print(f"WARNING: [{sample}] missing {label}: {candidate}", file=sys.stderr)
    return ""


def main():
    args = parse_args()

    if not args.mettanotator_results.is_dir():
        sys.exit(f"ERROR: mettanotator results directory not found: {args.mettanotator_results}")
    if args.virify_results is not None and not args.virify_results.is_dir():
        sys.exit(f"ERROR: virify results directory not found: {args.virify_results}")

    samples = read_input_samplesheet(args.input_samplesheet)

    # Validate that sample names actually match the mettanotator output structure.
    # A bare directory named after the sample is NOT enough: a wrong path can
    # coincidentally contain such a folder. Require that the sample directory
    # actually holds at least one expected mettanotator output file.
    matched = []
    for sample, assembly in samples:
        sample_dir = args.mettanotator_results / sample
        if not sample_dir.is_dir():
            print(
                f"WARNING: sample '{sample}' has no directory under "
                f"{args.mettanotator_results}; skipping it.",
                file=sys.stderr,
            )
            continue
        has_output = any(
            (sample_dir / rel.format(s=sample)).is_file()
            for rel in METTANOTATOR_PATHS.values()
        )
        if not has_output:
            print(
                f"WARNING: sample '{sample}' directory exists but contains no "
                f"expected mettanotator output files; skipping it.",
                file=sys.stderr,
            )
            continue
        matched.append((sample, assembly))

    # If nothing matched, the user almost certainly passed the wrong path(s).
    if not matched:
        sys.exit(
            "ERROR: none of the samples from the input samplesheet were found with "
            f"mettanotator output under {args.mettanotator_results}.\n"
            "       Check that --mettanotator-results points at the directory whose "
            "per-sample subfolders are named like the samplesheet's first column "
            "and contain the mettanotator annotation files."
        )

    samplesheet_rows = []
    manifest_rows = []
    n_virify = 0

    for sample, assembly in matched:
        mett = args.mettanotator_results

        proteins_gff = resolve(mett, METTANOTATOR_PATHS["proteins_gff"], sample, args.warn_missing, "proteins_gff (merged_gff)")
        proteins_faa = resolve(mett, METTANOTATOR_PATHS["proteins_faa"], sample, args.warn_missing, "proteins_faa (prokka)")
        ips_tsv      = resolve(mett, METTANOTATOR_PATHS["interproscan_tsv"], sample, args.warn_missing, "interproscan_tsv")
        amrfinder    = resolve(mett, METTANOTATOR_PATHS["amrfinder_tsv"], sample, args.warn_missing, "amrfinder_tsv")
        antismash    = resolve(mett, METTANOTATOR_PATHS["antismash_gff"], sample, args.warn_missing, "antismash_gff")
        gecco        = resolve(mett, METTANOTATOR_PATHS["gecco_gff"], sample, args.warn_missing, "gecco_gff")
        sanntis      = resolve(mett, METTANOTATOR_PATHS["sanntis_gff"], sample, args.warn_missing, "sanntis_gff")

        # virify gff is optional; empty when the dir isn't given or the file is absent
        virify_gff = ""
        if args.virify_results is not None:
            candidate = args.virify_results / sample / VIRIFY_REL.format(s=sample)
            if candidate.is_file():
                virify_gff = str(candidate.resolve())
                n_virify += 1

        samplesheet_rows.append({
            "sample": sample,
            "assembly": assembly,
            "proteins_gff": proteins_gff,
            "proteins_faa": proteins_faa,
            "virify_gff": virify_gff,
            "interproscan_tsv": ips_tsv,
        })
        manifest_rows.append({
            "sample": sample,
            "amrfinder_tsv": amrfinder,
            "antismash_gff": antismash,
            "gecco_gff": gecco,
            "sanntis_gff": sanntis,
        })

    # A virify directory was given but nothing resolved under it: likely a wrong
    # path. Warn (don't fail) since virify is optional and may legitimately have
    # no predictions for any sample.
    if args.virify_results is not None and n_virify == 0:
        print(
            f"WARNING: a virify results directory was given ({args.virify_results}) "
            "but no virify GFF was found for any sample. Check the path if you "
            "expected viral predictions.",
            file=sys.stderr,
        )

    _write_csv(args.out_samplesheet, MAP_SAMPLESHEET_COLS, samplesheet_rows)
    _write_csv(args.out_manifest, MANIFEST_COLS, manifest_rows)

    # ---- stdout summary ----
    n_input = len(samples)
    n = len(matched)
    print("MAP input generation complete.")
    print(f"  Samples in input samplesheet:   {n_input}")
    print(f"  Found in mettanotator results:  {n}")
    if args.virify_results is not None:
        print(f"  Found in virify results:        {n_virify}")
    else:
        print(f"  Found in virify results:        n/a (no virify path given)")
    print(f"  Samples written:                {n}")
    print(f"  With virify results:            {n_virify}")
    print(f"  Without virify results:         {n - n_virify}")
    print(f"  MAP samplesheet:                {args.out_samplesheet}")
    print(f"  Annotation manifest:            {args.out_manifest}")


def _write_csv(path, columns, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
