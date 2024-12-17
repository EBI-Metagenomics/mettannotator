process PSEUDOFINDER_GETDB {

    tag "Pseudofinder DB 2024_06"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: "copy"

    output:
    tuple path("uniprot_sprot/uniprot_sprot.fasta"), env("VERSION"), emit: pseudofinder_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pseudofinder/uniprot_sprot_2024_06.tar.gz

    mkdir -p uniprot_sprot

    tar -zxvf uniprot_sprot_2024_06.tar.gz -C uniprot_sprot

    VERSION=\$(cat uniprot_sprot/VERSION.txt)
    """
}
