
process AMRFINDER_PLUS_GETDB {

    tag "AMR Finder DB 2023-11-15.1"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.11.4--h6e70893_0'

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("amrfinder", type: "dir"), val("2023-11-15.1"), emit: amrfinder_plus_db

    script:
    """
    wget -r -nH --cut-dirs=5 \\
    ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/2023-11-15.1/

    amrfinder_index 2023-11-15.1

    mv 2023-11-15.1 amrfinder

    echo "2023-11-15.1" > amrfinder/VERSION.txt
    """
}
