
process BAKTA_GETDB {

    tag "BAKTA DB 2024-01-19"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("bakta", type: "dir"), val("2024-01-19"), emit: bakta_db

    script:
    """
    wget https://zenodo.org/record/10522951/files/db.tar.gz
    tar -xzf db.tar.gz
    rm db.tar.gz
    mv db bakta
    wget -r -nH --cut-dirs=5 \\
    ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/2023-11-15.1/
    rm -r bakta/amrfinderplus-db/*
    mv 2023-11-15.1 bakta/amrfinderplus-db/

    amrfinder_index bakta/amrfinderplus-db/2023-11-15.1

    echo "2024-01-19" > bakta/Version.txt
    """
}
