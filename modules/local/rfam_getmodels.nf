
process RFAM_GETMODELS {

    tag "Rfam models - release 14.9"

    container "quay.io/biocontainers/gnu-wget:1.18--h36e9172_9"

    publishDir "$params.dbs/rfam_models", mode: 'copy'

    output:
    tuple path("rfam_ncrna_cms"), val("14.9"), emit: rfam_ncrna_cms

    script:
    """
    mkdir -p rfam_ncrna_cms/

    wget --recursive --cut-dirs=5 -nH ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/ncrna/ -P rfam_ncrna_cms

    echo '14.9' > rfam_ncrna_cms/VERSION.txt
    """
}
