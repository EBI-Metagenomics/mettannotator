
process RFAM_GETMODELS {

    tag "Rfam models - release 15.0"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "$params.dbs/rfam_models", mode: 'copy'

    output:
    tuple path("rfam_ncrna_cms"), val("15.0"), emit: rfam_ncrna_cms

    script:
    """
    mkdir -p rfam_ncrna_cms/

    wget --recursive --cut-dirs=5 -nH ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/rfam_15.0/ncrna_cms/ -P rfam_ncrna_cms

    echo '15.0' > rfam_ncrna_cms/VERSION.txt
    """
}
