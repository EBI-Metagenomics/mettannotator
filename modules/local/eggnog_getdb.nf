process EGGNOG_MAPPER_GETDB {

    tag "EGGNOG Mapper DB 5.0.2"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9' :
        'biocontainers/gnu-wget:1.18--h36e9172_9' }"

    publishDir "${params.dbs}", mode: "copy"

    output:
    tuple path("eggnog/", type: "dir"), val("5.0.2"), emit: eggnog_db

    script:
    """
    mkdir -p eggnog/

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/eggnog_db_5.0.2.tgz

    tar -xvzf eggnog_db_5.0.2.tgz -C eggnog/

    echo "5.0.2" > eggnog/VERSION.txt

    rm eggnog_db_5.0.2.tgz
    """
}
