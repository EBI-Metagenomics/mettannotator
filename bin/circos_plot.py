#!/usr/bin/env python3
# coding=utf-8

import argparse
import logging
import sys

from pycirclize import Circos
from pycirclize.parser import Gff

from matplotlib.patches import Patch

logging.basicConfig(level=logging.INFO)


def main(infile, outfile, prefix, mobilome):
    modified_infile = remove_escaped_characters(infile)
    gff = Gff(modified_infile)
    seqid2size = gff.get_seqid2size()
    if len(seqid2size) > 200:
        logging.info(
            "Skipping plot generation for file {} due to a large number of contigs: {}. "
            "Plots are only generated for genomes with up to 200 annotated contigs.".format(
                infile, len(seqid2size)
            )
        )
        sys.exit()

    seqid2features = gff.get_seqid2features(feature_type=None)

    circos = Circos(seqid2size, space=2, start=1, end=358)

    circos.text("{}\n".format(prefix), size=15, r=30)

    for sector in circos.sectors:
        # Plot contig labels
        sector.text(
            sector.name, orientation="vertical", r=110, size=10, color="dimgrey"
        )
        # Plot scale
        position_track = sector.add_track((99, 100))
        position_track.axis(fc="none")
        major_ticks_interval = 500000
        minor_ticks_interval = 50000
        if sector.size > minor_ticks_interval:
            position_track.xticks_by_interval(
                major_ticks_interval, label_formatter=lambda v: f"{v / 10 ** 6:.1f} Mb"
            )
            position_track.xticks_by_interval(
                minor_ticks_interval, tick_length=1, show_label=False
            )

        # Initiate feature tracks
        f_cds_track = sector.add_track((93, 98), r_pad_ratio=0.1)
        f_cds_track.axis(fc="none", ec="none")
        r_cds_track = sector.add_track((88, 93), r_pad_ratio=0.1)
        r_cds_track.axis(fc="none", ec="none")
        rna_track = sector.add_track((83, 87), r_pad_ratio=0.1)
        rna_track.axis(fc="none", ec="none")
        bgc_track_antismash = sector.add_track((78, 80), r_pad_ratio=0.1)
        bgc_track_antismash.axis(fc="none", ec="tomato", ls="dashdot", lw=0.15)
        bgc_track_gecco = sector.add_track((76, 78), r_pad_ratio=0.1)
        bgc_track_gecco.axis(fc="none", ec="lightsalmon", ls="dashdot", lw=0.15)
        bgc_track_sanntis = sector.add_track((74, 76), r_pad_ratio=0.1)
        bgc_track_sanntis.axis(fc="none", ec="firebrick", ls="dashdot", lw=0.15)
        dbcan_track = sector.add_track((68, 70), r_pad_ratio=0.1)
        dbcan_track.axis(fc="none", ec="forestgreen", ls="dashdot", lw=0.15)
        amr_track = sector.add_track((62, 64), r_pad_ratio=0.1)
        amr_track.axis(fc="none", ec="dodgerblue", ls="dashdot", lw=0.15)
        antiphage_track = sector.add_track((56, 58), r_pad_ratio=0.1)
        antiphage_track.axis(fc="none", ec="orchid", ls="dashdot", lw=0.15)
        if mobilome:
            mobilome_track = sector.add_track((50, 52), r_pad_ratio=0.1)
            mobilome_track.axis(fc="none", ec="lightseagreen", ls="dashdot", lw=0.15)

        for feature in seqid2features[sector.name]:
            if feature.type == "CDS":
                if feature.strand == 1:
                    f_cds_track.genomic_features([feature], fc="hotpink")
                else:
                    r_cds_track.genomic_features([feature], fc="steelblue")
                if "antismash_bgc_function" in feature.qualifiers:
                    bgc_track_antismash.genomic_features([feature], fc="tomato")
                if "gecco_bgc_type" in feature.qualifiers:
                    bgc_track_gecco.genomic_features([feature], fc="lightsalmon")
                if "nearest_MiBIG" in feature.qualifiers:
                    bgc_track_sanntis.genomic_features([feature], fc="firebrick")
                if "dbcan_prot_type" in feature.qualifiers:
                    dbcan_track.genomic_features([feature], fc="forestgreen")
                if "amrfinderplus_scope" in feature.qualifiers:
                    amr_track.genomic_features([feature], fc="dodgerblue")
                if "defense_finder_type" in feature.qualifiers:
                    antiphage_track.genomic_features([feature], fc="orchid")

            elif feature.type in ["tRNA", "ncRNA"]:
                rna_track.genomic_features([feature], fc="darkmagenta")
            elif mobilome and feature.type in [
                "mobility_island",
                "cellular_recombinase",
                "insertion_sequence",
                "conjugative_element",
                "nested_mobile_element",
                "terminal_inverted_repeat_element",
                "direct_repeat_element",
            ]:
                mobilome_track.genomic_features([feature], fc="lightseagreen")
            elif mobilome and feature.type.lower() == "phage":
                mobilome_track.genomic_features([feature], fc="blue")

    fig = circos.plotfig()
    # Add legend
    handles = [
        Patch(color="hotpink", label="Forward CDS"),
        Patch(color="steelblue", label="Reverse CDS"),
        Patch(color="darkmagenta", label="RNA"),
        Patch(color="tomato", label="BGCs (antiSMASH)"),
        Patch(color="lightsalmon", label="BGCs (GECCO)"),
        Patch(color="firebrick", label="BGCs (SanntiS)"),
        Patch(color="forestgreen", label="Predicted PULs"),
        Patch(color="dodgerblue", label="AMR genes"),
        Patch(color="orchid", label="Anti-phage defense genes"),
    ]
    if mobilome:
        handles = handles + [
            Patch(color="blue", label="Mobilome (phage)"),
            Patch(color="lightseagreen", label="Mobilome (other)"),
        ]

    main_legend = circos.ax.legend(
        handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8
    )
    circos.ax.add_artist(main_legend)
    fig.savefig(outfile, dpi=600)


def remove_escaped_characters(infile):
    outfile = infile + "_modified"
    with open(infile, "r") as file_in:
        content = file_in.read()
        modified_content = content.replace("\\=", "")

    # Open the file in write mode to overwrite it with the modified content
    with open(outfile, "w") as file_out:
        file_out.write(modified_content)
    return outfile


def parse_args():
    parser = argparse.ArgumentParser(description="Test script for Circos plots.")
    parser.add_argument(
        "-i", "--infile", required=True, help="Path to the GFF file to plot"
    )
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the output file"
    )
    parser.add_argument(
        "-p", "--prefix", required=True, help="Prefix to use for the genome"
    )
    parser.add_argument(
        "--mobilome",
        required=False,
        action="store_true",
        default=False,
        help="Plot the mobilome track. Default: False",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args.infile, args.outfile, args.prefix, args.mobilome)
