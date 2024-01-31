
process RFAM_GETMODELS {

    tag "RFam Models"

    container "quay.io/biocontainers/gnu-wget:1.18--h36e9172_9"

    output:
    path "rfams_cms", emit: rfams_cms
    path "ncrna_cms", emit: ncrna_cms

    script:
    """
    mkdir rfams_cms/

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/rfams_cms/ -P rfams_cms

    mkdir ncrna_cms/

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/ncrna/ -P ncrna_cms
    """
}
