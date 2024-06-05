
process BAKTA_GETDB {

    tag "BAKTA DB 2024-01-19"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'

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
    ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/
    rm -r bakta/amrfinderplus-db/*
    mv 2024-01-31.1 bakta/amrfinderplus-db/

    amrfinder_index bakta/amrfinderplus-db/2024-01-31.1

    mv bakta/amrfinderplus-db/2024-01-31.1/ bakta/amrfinderplus-db/latest/

    echo "2024-01-19" > bakta/VERSION.txt
    """
}
