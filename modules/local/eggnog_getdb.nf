
process EGGNOG_MAPPER_GETDB {

    tag "EGGNOG Mapper DB 5.0.2"

    container "quay.io/biocontainers/gnu-wget:1.18--h36e9172_9"

    output:
    path "eggnog_data/eggnog.db"           , emit: eggnog_db
    path "eggnog_data/eggnog_proteins.dmnd", emit: eggnong_diamond_db
    path "eggnog_data/", type: "dir"       , emit: eggnong_data_dir

    script:
    """
    mkdir eggnog_data/

    wget ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/genomes-pipeline/eggnog_db_5.0.2.tgz

    tar -xvzf eggnog_db_5.0.2.tgz -C eggnog_data/

    rm eggnog_db_5.0.2.tgz
    """
}
