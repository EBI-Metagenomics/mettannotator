process PSEUDOFINDER_GETDB {

    tag "Pseudofinder DB 2024_06"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: "copy"

    output:
    tuple path("uniprot_sprot_2024_06.fasta"), val("2024_06"), emit: pseudofinder_db

    script:
    """
    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/pseudofinder/uniprot_sprot_2024_06.fasta.gz

    gunzip uniprot_sprot_2024_06.fasta.gz

    """
}
