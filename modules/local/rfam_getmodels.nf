
process RFAM_GETMODELS {

    tag "RFam Models"

    container "quay.io/biocontainers/gnu-wget:1.18--h36e9172_9"

    publishDir "$params.dbs/rfam_models", mode: 'copy'

    output:
    path "rfam_rrna_cms", emit: rfam_rrna_cms
    path "rfam_ncrna_cms", emit: rfam_ncrna_cms

    script:
    """
    mkdir -p rfam_rrna_cms/

    wget --recursive --cut-dirs=5 -nH ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/rfams_cms/ -P rfam_rrna_cms

    mkdir -p rfam_ncrna_cms/

    wget --recursive --cut-dirs=5 -nH ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/ncrna/ -P rfam_ncrna_cms
    """
}
