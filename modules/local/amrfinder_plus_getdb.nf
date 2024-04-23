
process AMRFINDER_PLUS_GETDB {

    tag "AMR Finder DB 2024-01-31.1"

    container 'quay.io/biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0'

    publishDir "$params.dbs", mode: 'copy'

    output:
    tuple path("amrfinder", type: "dir"), val("2024-01-31.1"), emit: amrfinder_plus_db

    script:
    """
    wget -r -nH --cut-dirs=5 \\
    ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.12/2024-01-31.1/

    amrfinder_index 2024-01-31.1

    mv 2024-01-31.1 amrfinder

    echo "2024-01-31.1" > amrfinder/VERSION.txt
    """
}
