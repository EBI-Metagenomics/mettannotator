#!/usr/bin/env python3
# coding=utf-8

import argparse
import logging
import sys

from pycirclize import Circos
from pycirclize.parser import Gff

from matplotlib.patches import Patch

logging.basicConfig(level=logging.INFO)


def main(infile, outfile, prefix, contig_num_limit, contig_trim, mobilome, dpi):
    modified_infile = remove_escaped_characters(infile)
    gff = Gff(modified_infile)
    seqid2size = gff.get_seqid2size()
    if len(seqid2size) > contig_num_limit:
        logging.info(
            "Skipping plot generation for file {} due to a large number of contigs: {}. "
            "Plots are only generated for genomes with up to {} annotated contigs.".format(
                infile, len(seqid2size), contig_num_limit
            )
        )
        sys.exit()

    seqid2features = gff.get_seqid2features(feature_type=None)

    circos = Circos(seqid2size, space=1, start=1, end=358)

    circos.text("{}\n".format(prefix), size=15, r=30)

    # Skip printing contig names if the names are too long
    print_contigs = True
    if contig_trim == 500:  # the user didn't choose truncation, check lengths
        for sector in circos.sectors:
            if len(sector.name[:contig_trim]) > 24:
                print_contigs = False
    if not print_contigs:
        logging.info("Not printing contig labels because they are too long. Rerun the script with the "
                     "--contig-trim flag to truncate the labels if you would like them printed.")
    for sector in circos.sectors:
        if print_contigs:
            # Plot contig labels
            sector.text(
                sector.name[:contig_trim], orientation="vertical", r=110, size=6, color="dimgrey"
            )
        # Plot scale
        position_track = sector.add_track((99, 100))
        position_track.axis(fc="none")
        major_ticks_interval = 500000
        minor_ticks_interval = 50000
        if sector.size > minor_ticks_interval:
            if sector.size >= minor_ticks_interval * 10:
                position_track.xticks_by_interval(
                    major_ticks_interval, label_formatter=lambda v: f"{v / 10 ** 6:.1f} Mb"
                )
            else:
                position_track.xticks_by_interval(
                    major_ticks_interval, show_label=False
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
                "viral_sequence",
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
    fig.savefig(outfile, dpi=dpi)


def remove_escaped_characters(infile):
    outfile = infile + "_modified"
    with open(infile, "r") as file_in:
        content = file_in.read()
        modified_content = content.replace("\\=", "")

    # Write the modified content into a file
    with open(outfile, "w") as file_out:
        file_out.write(modified_content)
    return outfile


def parse_args():
    parser = argparse.ArgumentParser(description="Script for Circos plot generation.")
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
        "-l",
        "--limit",
        required=False,
        type=int,
        default=50,
        help="Only generate a plot if the genome has no more than this number of contigs. Limit introduced because "
             "highly fragmented genomes do not produce readable plots. Default: 50.",
    )
    parser.add_argument(
        "--contig-trim",
        required=False,
        default=500,
        type=int,
        help="If the contig length is over 24 characters long, contig names will not be printed on the plot. Specify "
             "the length to trim the contig names down to if you would like the shorter names printed.",
    )
    parser.add_argument(
        "--mobilome",
        required=False,
        action="store_true",
        default=False,
        help="Plot the mobilome track. Default: False",
    )
    parser.add_argument(
        "--dpi",
        required=False,
        type=int,
        default=600,
        help="Specify the dpi for the plot. Default: 600",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args.infile, args.outfile, args.prefix, args.limit, args.contig_trim, args.mobilome, args.dpi)
